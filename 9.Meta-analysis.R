# sn-RNA-seq sex-difs atlas in PD
# Meta-analysis
# Autor: Fernando Gordillo Gonzalez
# Args: see below

## starting
date ()
Sys.info ()[c("nodename", "user")]
R.version.string
setwd("~/PDsex_atlas")

#-------------------------------- Arguments --------------------------------####
library(argparse)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-d", "--dir", type="character", default="/results/8.Differential_exp/MAST/", 
                    help="Directory in each study folder that contains DEA results")
parser$add_argument("-c", "--contrast", type="character", default="FM", 
                    help="MA between case and control(PDvsC), female(F), male (M) or sex differences(FM)")
parser$add_argument("-ct", "--cell", type="character", default=NULL, 
                    help="Cell type to meta-analyze")
parser$add_argument("-t", "--tissue", type="character", default="SN", 
                    help="Tissue selected")
parser$add_argument("-m", "--metadata", type="character", default="Integracion/studies.tsv", 
                    help="File containing studies tissue (tsv)")
parser$add_argument("-p", "--pvalue", type="double", default=0.05, 
                    help="Pvalue selected as threshold for the analysis, defaults at 0.05")
parser$add_argument("-o", "--outputdir", type="character", default="Integracion/MA/", 
                    help="Output directory to save results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

if (!file.exists(args$metadata) | strsplit(args$metadata,"\\.")[[1]][2] != "tsv"){
  stop(paste(args$metadata,"not exist or not supported (tsv)"))
} else { 
  metastudies=read.delim(args$metadata)
  stopifnot("Not suported contrast."=args$contrast %in% c("PDvsC","F", "M","FM"),
            "Tissue not among possibles (Tissue column from metadata file)."=args$tissue %in% c(metastudies$Tissue,"all"),
            "Not all directories exists. Check --dir param or metadata file."={unlist(lapply(metastudies$Studies, FUN=function(x){dir.exists(paste0(x,args$dir))}))})
}

#----------------------------- Load packages -------------------------------####

library(Biobase)
library(metafor)
library(tibble)
library(limma)
library(dplyr)

#---------------------- Make directories & Load data -----------------------####
metaanalysis_name <- paste("MA", args$tissue,args$cell ,args$contrast, sep="_")
print(metaanalysis_name)
meta_dir <- paste0(getwd(),"/",args$outputdir ,args$tissue,"/",args$cell ,"/",args$contrast)
dir.create(meta_dir, recursive = TRUE)

## load data
selected_studies=if(args$tissue == "all"){ metastudies$Studies 
                 } else { metastudies$Studies[metastudies$Tissue==args$tissue] }
       
EDs=list()
for (i in selected_studies){
  print(paste("Loading study", i))
  EDs[[i]]=read.delim(paste0(i,args$dir,args$cell ,"/all_res",args$contrast,".tsv"), row.names=1)
}

setwd(meta_dir)


#---------------------------- Meta-analysis --------------------------------####

# STEP 1. Preparing input for meta-analysis: LOR and SE matrix ####
# ===============================================================

# we search a list including all unique ID genes for all studies
genes <- unlist(lapply(EDs, rownames))
genes <- unique (genes)
genes <- sort (genes)


### generating matrix with all SE for all studies
mat.SE <- matrix (NA, nrow = length (genes), ncol = length(EDs))
rownames (mat.SE) <- genes
colnames (mat.SE) <- names(EDs)

for (i in 1:length(EDs)){
  mat.SE[, names(EDs[i])] <- EDs[[i]][rownames(mat.SE), "SE"]
}
mat.SE[mat.SE==0]=NA
cat("\nInitial dims matrix SE:")
dim (mat.SE)

# select genes included at least in 2 or more studies
mat.SE.NA <- is.na(mat.SE)
sum.NA <-  apply(mat.SE.NA, 1, sum)
min.sum.NA <- sum.NA < ncol(mat.SE) - 1

# filter by min.sum.NA
mat.SE <- mat.SE[min.sum.NA == T, ]
mat.SE.NA <- mat.SE.NA[min.sum.NA == T, ]

# filter feats that does not appear in both tissues studied
# if(args$tissue == "all"){
#   not_all <- rep(FALSE, nrow(mat.SE.NA))
#   for (t in metastudies$Tissue){
#     tistudies=metastudies$Studies[metastudies$Tissue==t]
#     not_t=apply(mat.SE.NA[,tistudies], 1,sum)==length(tistudies)
#     not_all <- not_all | not_t
#   }
#   mat.SE <- mat.SE[!not_all, ]
# }

cat("\nFinal dims matrix SE:")
dim (mat.SE)

### generating matrix with all logFC for all studies
mat.logFC <- matrix (NA, nrow = length(genes), ncol = length(EDs))
rownames (mat.logFC) <- genes
colnames (mat.logFC) <- names(EDs)

for (i in 1:length(EDs)){
  mat.logFC[, names(EDs[i])] <- EDs[[i]][rownames(mat.logFC), "logFC"]
}

cat("\nInitial dims matrix logFC:")
dim (mat.logFC)

# filter by min.sum.NA
mat.logFC <- mat.logFC[min.sum.NA == T, ]

# filter feats that does not appear in both tissues studied
# if(args$tissue == "all"){
#   mat.logFC <- mat.logFC[!not_all, ]
# }
cat("\nFinal dims matrix logFC:")
dim (mat.logFC)

# STEP 2. Meta-analysis for genes ####
# ===============================================================

# suppose between-study variance is non-zero.
# there are different methods to estimate this variance:
# DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
# Now we have logFC and SE  (not VARIANCE), so:
# yi -> logFC   sei -> SE
# result.lor <- rma(yi = mat.logFC[1, ],
#                   sei = mat.SE[1, ],   #pay attention, not vi (varianze)
#                   method = "DL") # DerSimonian-Laird.

# explore the function to do the meta-analysis
#?rma

MA <- lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SE[x, ],
                             method = "DL")})

names (MA) <- rownames(mat.logFC)

print(paste0("Dim meta-analysis: ", length(MA)))

#data.frame including all detailed results:
result_meta <- as.data.frame(do.call("rbind",
                                     lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub,
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                            })))

colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                           "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")

result_meta$p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr")
result_meta$p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
result_meta$logFDR <- -log10(result_meta$p.adjust.fdr)

# significant genes
# En el caso de SexDiffs no sale significativos con .05, lo triplicamos
print(paste("Genes significant (p<",args$pvalue ,"):",sum(result_meta[, "pvalue"] < args$pvalue )))
print(paste("Genes significant (FDR<",args$pvalue ,"):",sum(result_meta[, "p.adjust.fdr"] < args$pvalue )))

# add number of studies where the gene is evaluated
n.studies <-  ncol(mat.logFC) - sum.NA
table(n.studies)
result_meta$n.studies <- n.studies[rownames(mat.logFC)]

sig.genes.df = result_meta[result_meta$p.adjust.fdr < args$pvalue ,]

write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),] %>% rownames_to_column('ENSEMBL'), file = "sig.genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = result_meta[order(result_meta$p.adjust.fdr),] %>% rownames_to_column('ENSEMBL'), file = "all.genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS ####
# ===============================================================

#add 4 new variables about influence & sensitivity analysis:

for (i in rownames(sig.genes.df)){
  #print(i)
  #define studies for each function (not NA)
  estudios <- colnames(mat.logFC)[!mat.SE.NA[i,]]

  #influence info 1:
  #number of studies where the sign of the logOR is the same  of the global logOR:
  sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios)))

  #influence info 2: how many studies could be influencers?
  inf <- influence(MA[[i]])
  res <- paste(estudios[inf$is.infl], collapse = ",")
  sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)

  #sensivity analysis
  l1 <-as.data.frame(leave1out(MA[[i]]))
  rownames(l1) <- estudios

  #1. p.value about differences between all estimates from leave one out
  #   and global estimate)
  sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                           mu=as.numeric(MA[[i]]$b))$p.value
  #2. number of  studies where pvalue > 0.05
  # (we hope p-values < 0.05, significant estimates)
  res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")
  sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}

## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:

#1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
cat("\n\n Individual LORs with same sign as global LOR \n")
table(sig.genes.df$infl.same.sign.logFC)

#2. INFLUENCE STUDIES. How many functions including influence studies?
cat("\n\n Genes with influential studies \n")
table(sig.genes.df$infl.nstudies=="non")

#3. SENSITIVITY. In global, are there many functions with differences in the estimate?
cat("\n\n leave 1 out (sensi.global < 0.05) \n")
table(sig.genes.df$sensi.global < 0.05)

#4. SENSITIVITY.  How many functions including changes in the significance about
# its new estimate  after leave1out?
cat("\n\n leave 1 out doesnt change p value \n")
table(sig.genes.df$sensi.specific == "all.p.values < 0.05")


#save final results:

write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),] %>% rownames_to_column('ENSEMBL'), 
            file = "sig.genes_def.tsv", sep ="\t", quote = F, row.names = F)


# STEP 4. Visualization of significant genes ####
# ===============================================================
sig.genes.df.5 <- head(sig.genes.df[order(sig.genes.df$p.adjust.fdr),], n=5)
#select significant functions to visualize:
sig.results <- result_meta[rownames(sig.genes.df.5),]

dir.create("./plots")
setwd("plots")

selMethod <- "DL"

for (i in 1:nrow(sig.results)){
  mygenes <- rownames(sig.results)[i]
  res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")

  #FOREST PLOT
  # png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 ,
  #      height = 960, res = 200)
  png (filename = paste(gsub("-","_",mygenes),"_FOREST",".png", sep =""), width = 960 ,
       height = 960, res = 200)
  forest(res,
         slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
         xlab="logFC", cex=0.7,
         mlab=paste(selMethod, "Model for All Studies", sep = " "),
         border = "black", #Color del borde del rombo
         col = "red", #Color del rombo
         main = paste("\n", mygenes, sep=""))
  text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
  dev.off()

  #FUNNEL PLOT
  png (filename = paste(gsub("-","_",mygenes),"_FUNNEL", ".png", sep =""), width = 960 ,
       height = 960, res = 200)
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance",
         back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
  par(mfrow=c(1,1))
  dev.off()

  #INFLUENCE PLOTS
  # That shows various diagnostic measures
  png (filename = paste(gsub("-","_",mygenes), "_INFLUENCE", ".png", sep =""), width = 960 ,
       height = 960, res = 200) ##CAMBIAR
  inf <- influence(res)
  #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
  plot(inf)
  dev.off()

}
print("Generating report")

# STEP 5. Generating report ####
# ===============================================================


# Summary table
res <- function(df){
  total <- nrow(df)
  sig <- sum(df$p.adjusted < 0.05)
  res <- c(total, sig)
  names(res) <- c("genes", "sig.genes")
  return(res)
}

tabla_resumen <- function(x) {
  table <- lapply(x, res)
  # Las uno en una tabla
  table <- dplyr::bind_rows(table, .id = "Study")

    return(table)
}

#Tabla resumen con los genes y los genes significativos (p.val.adj < 0.05)
table_sum <- tabla_resumen(EDs)
# table_sum$Study <- gsub("_ED","", table_sum$Study)
# if(all(names(apply(!is.na(mat.logFC), 2, sum)) == table_sum$Study)){
#   table_sum <- dplyr::bind_cols(table_sum, apply(!is.na(mat.logFC), 2, sum))
#   table_sum <- table_sum %>% dplyr::relocate(...4, .after = genes)
# }

# Function to create multiple tabs
make.tabs <- function(sig.genes.df){
  res <- NULL
  for(g in rownames(sig.genes.df)){
    file_name <- gsub("-","_", g)
    res <- c(res, '## ', g, ' {.tabset .tabset-fade .tabset-pills} \n',
             "**Statistics of ", g, " meta-analisys** \n",
             "```{r, fig.align='center'}", '\n',
             "kable(sig.genes.df['",g,"',])", '\n',
             '```', '\n\n',
             "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sig.genes.df[g,"SYMBOL"], ") \n\n",
             "### Forest plot {-} \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_FOREST.png")\n',
             '```', ' \n\n',
             "### Funnel plot {-} \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_FUNNEL.png")\n',
             '```', ' \n\n',
             "### Incluence plot {-} \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_INFLUENCE.png")\n',
             '```', '\n\n')
  }
  return(res)
}


# Generating volcano plotly
sig_limit <- 0.05
lfc_limit <- 0.5

data <- result_meta
data$gene_name <- data$SYMBOL
rownames(data) <- NULL
data <- data[,c("gene_name", "logFC", "p.adjust.fdr", "logFDR")]

# add a grouping column; default value is "not significant"
data["Group"] <- "NotSignificant"

# change the grouping for the entries with significance but not a large enough Fold change
data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) < lfc_limit ),"Group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
data[which(data['p.adjust.fdr'] > sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("FoldChange > ",lfc_limit," | < -",lfc_limit)

# change the grouping for the entries with both significance and large enough fold change
data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("Significant & FoldChange > ",lfc_limit," | < -",lfc_limit)

data$texto <- paste0("Gene: ", data$gene_name, "\np.val.fdr: ", data$p.adjust.fdr)

library(ggplot2)
volcano_plot <- ggplot(data, aes(x = logFC, y = logFDR)) +
  geom_point(aes(color = Group, text = texto)) +
  scale_color_manual(values = c("#E69F00", "#DB5E00", "#0072B2", "#009E73")) + #https://stackoverflow.com/q/57153428
  xlab("logFC") + ylab("-log10(FDR)")

# Create the Rmd to knit
ns <- nrow(sig.genes.df)
cat(
  '---
title: "Global Meta-analysis"
output:
  html_document:
    toc: false
    toc_float: false
    code_folding: hide
    number_sections: true
    theme: spacelab
---
## ', metaanalysis_name, ' {.tabset .tabset-pills -}

```{r, warning=F, message=F}
library(dplyr)
library(knitr)
library(DT)
library(ggplot2)
library(plotly)
load("', meta_dir, '/', metaanalysis_name, '.RData")
```  \n
### MA results (significant genes) {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "datatable(sig.genes.df[,c(1:5,12,15)], caption='Statistics of ", ns, " significant genes',rownames=FALSE, escape = FALSE, filter = 'top')", '\n',
  '```', '\n\n',
  '### MA results (volcano plot) {-} \n',
  "```{r, fig.align='center', out.width = '99%'}", '\n',
  "ggplotly(volcano_plot)", '\n',
  '```', '\n\n',
  '### Individual results (summary) {-} \n',
  "```{r, fig.align='center'}", '\n',
  "kable(table_sum, col.names = c('Study','All genes','Significant genes'))", '\n',
  '```',
  '\n\n',
  '### sessionInfo {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "date()", "\n",
  "sessionInfo()", '\n',
  '```', '\n\n',
  sep = "",
  file = paste0(meta_dir, "/", metaanalysis_name, ".Rmd"))

cat(
'---
title: "Meta-analysis of genes"
subtitle: "Draft, ', format(Sys.time(), "%d-%m-%Y"),'v6c"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    code_folding: hide
    number_sections: false
    theme: spacelab
---
# ', metaanalysis_name, ' {- .tabset .tabset-fade .tabset-pills}

```{r, warning=F, message=F}
library(dplyr)
library(knitr)
library(DT)
load("', meta_dir, '/', metaanalysis_name, '.RData")
```  \n',
make.tabs(sig.genes.df.5[order(row.names(sig.genes.df.5)),]), "\n\n",
'## sessionInfo {-}  \n',
"```{r, fig.align='center'}", '\n',
"date()", "\n",
"sessionInfo()", '\n',
'```', '\n\n',
  sep = "",
  file = paste0(meta_dir, "/", metaanalysis_name, "_genes.Rmd"))

save(sig.genes.df, result_meta, table_sum, volcano_plot, file = paste0(meta_dir, "/", metaanalysis_name, ".RData"))
# Render the Rmd created into html here ####
rmarkdown::render(paste0(meta_dir, "/", metaanalysis_name, ".Rmd"))
rmarkdown::render(paste0(meta_dir, "/", metaanalysis_name, "_genes.Rmd"))
