# sn-RNA-seq sex-difs atlas in PD
# Differential expression analysis (MAST)
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID cell_type gene_lengths_path

#----------------------------- Load packages ----------------------------------####

library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(MAST)
library(viridis)
library(BiocParallel)

#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
cell_type = as.character(args[2])

# Create directories
dir.create(paste0(estudio,"/results/8.Differential_exp/MAST/",cell_type), recursive = T)
setwd(paste0(estudio,"/results/8.Differential_exp/MAST/",cell_type))
getwd()

# Load data
sce = readRDS(paste0("../../../5.Cell_annotation/sce_",tolower(estudio),"_annot_filt.rds"))

# Selecting cell type
sce = sce[,sce$cell_type == cell_type]

#------------------------------- Functions ---------------------------------####

## TPMs function
tpm.local = function(counts, len) {
  x = counts/len
  return(t(t(x)*1e6/colSums(x)))
}

zz_process = function(zlm, sce) {
  
  results = cbind(data.frame(rownames(zlmCond_cells)), data.frame(zlmCond_cells[,3,c(1,3)]))
  colnames(results) = c("gene_id","lambda","p.value")
  results$symbol=rowData(sce)[match(results$gene_id, rownames(sce)), "SYMBOL"]
  results=results[,c(1,4,2,3)]
  results[,"p.adjusted"] = p.adjust(results$p.value, "BH") # Correccion del pvalor
  
  return(results)
}

logFC_mast = function(zlm, results, case, control) {
  coefnames = colnames(zlm@LMlike@modelMatrix)
  
  # Constrast 0
  contrast0 = setNames(rep(0, length(coefnames)), coefnames)
  contrast0[control] = 1
  contrast1=data.frame("case"=setNames(rep(0, length(coefnames)), coefnames))
  contrast1[case,] = 1
  
  # control=Hypothesis(control, colnames(zlm@LMlike@modelMatrix))
  # case=Hypothesis(case, colnames(zlm@LMlike@modelMatrix))
  
  fc=getLogFC(zlm, contrast0, contrast1)
  sum =summary(zlm, logFC=fc)$datatable
  sum=sum[contrast == "case" & component=='logFC',-c(2,3)]
  sum$SE <- (sum[, "ci.hi"] - sum[, "ci.lo"])/ 3.92
  colnames(sum)[c(1,4)]=c("gene_id","logFC")
  sum[is.na(sum)]=0
  
  results = merge(results, sum, by ="gene_id", sort=F)
  
  return(results)
}

logFC_cpm = function(sce, genes, case, control) {
  
  # Calculate cpms and select data of each group. Sum pseudocount factor (1)
  cpm_case = edgeR::cpm(counts(sce))[genes,sce$cond_sex==case]  + 1
  cpm_control = edgeR::cpm(counts(sce))[genes,sce$cond_sex==control]  + 1
  
  # Calculate logFC
  LFC = as.data.frame(log2(rowMeans(cpm_case)/rowMeans(cpm_control)))
  colnames(LFC) = "logFC.cpm"
  LFC$gene_id = rownames(cpm_case)
  
  return(LFC)
}

#------------------------------- Preprocesing ---------------------------------####

### Compute gene lengths
human.genes = read.delim(args[3], header = TRUE, sep = "\t", row.names = 1)
human.genes$gene.length=(human.genes$Gene.end..bp. - human.genes$Gene.start..bp.) + 1

idx<-match(rownames(sce),rownames(human.genes))
rowData(sce)$gene.length=human.genes[idx,"gene.length"]
idx_symb<-match(rowData(sce)[is.na(idx),]$SYMBOL,human.genes$Gene.name)
rowData(sce)[is.na(idx),]$gene.length=human.genes[idx_symb,"gene.length"]
sce=sce[!is.na(rowData(sce)$gene.length),]

### Compute log2(TPM + 1) values
tpm.counts = tpm.local(counts = assay(sce, "counts"), len = rowData(sce)[,"gene.length"])
assay(sce, "log2tpm") = log2(tpm.counts + 1)
rm(tpm.counts, idx, idx_symb)
# saveRDS(sce,paste0("sce_",tolower(estudio),"_annot_filt2.rds"))
#------------------------------ MAST analysis ---------------------------------####

print(paste0("Expresión diferencial con MAST en: ", cell_type))

# 1. Scaling number of genes
cdr2_cells = colSums(assay(sce, "log2tpm")>0)
colData(sce)$cngeneson = scale(cdr2_cells)

# 2. Hurdle model
# Construct sca object
rowData(sce)$primerid = rownames(sce)
sce.assay_cells = SceToSingleCellAssay(sce,check_sanity = FALSE)

# Model
zlmCond = zlm(~0+cond_sex+cngeneson+Sample, sce.assay_cells, exprs_values = "log2tpm", method = "bayesglm", ebayes = TRUE, parallel = TRUE)
saveRDS(zlmCond, "zlm_model.rds")

### CaseFemale - ControlFemale ####
# Test hypothesis
zlmCond_cells <- lrTest(zlmCond, Hypothesis("cond_sexPD_F - cond_sexC_F"))
saveRDS(zlmCond, "zlm_lrt_F.rds")

results = zz_process(zlmCond_cells,sce)
write.table(results, file=paste0("./all_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond,results,"cond_sexPD_F", "cond_sexC_F")
write.table(results, file=paste0("./all_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

logFC_F=logFC_cpm(sce, results[,"gene_id"] ,"PD_F", "C_F")
results = merge(results, logFC_F, by ="gene_id", sort=F)

results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

# Number of significative genes
sig = results[results$p.adjusted<0.05,]
write.table(sig, file=paste0("./sig_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)


### CaseMale - ControlMale ####
# Test hypothesis
zlmCond_cells <- lrTest(zlmCond, Hypothesis("cond_sexPD_M - cond_sexC_M"))
saveRDS(zlmCond, "zlm_lrt_M.rds")

results = zz_process(zlmCond_cells,sce)
write.table(results, file=paste0("./all_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond,results,"cond_sexPD_M", "cond_sexC_M")
write.table(results, file=paste0("./all_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

logFC_M=logFC_cpm(sce,results[,"gene_id"],"PD_M", "C_M")
results = merge(results, logFC_M, by ="gene_id", sort=F)

results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

# Number of significative genes
sig = results[results$p.adjusted<0.05,]
write.table(sig, file=paste0("./sig_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)


### (CaseFemale - ControlFemale) - (CaseMale - ControlMale) ####
# Test hypothesis
zlmCond_cells <- lrTest(zlmCond, Hypothesis("(cond_sexPD_F - cond_sexC_F) - (cond_sexPD_M - cond_sexC_M)"))
saveRDS(zlmCond, "zlm_lrt_FM.rds")

results = zz_process(zlmCond_cells,sce)
write.table(results, file=paste0("./all_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond, results, c("cond_sexC_M","cond_sexPD_F"), c("cond_sexC_F","cond_sexPD_M"))
write.table(results, file=paste0("./all_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

# logFC CPM
logFC_FM=merge(logFC_F, logFC_M, by ="gene_id")
logFC_FM$logFC.cpm = logFC_FM[,"logFC.cpm.x"] - logFC_FM[,"logFC.cpm.y"]
logFC_FM=logFC_FM[,c("gene_id","logFC.cpm")]
results = merge(results, logFC_FM, by ="gene_id", sort=F)


results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

# Number of significative genes
sig = results[results$p.adjusted<0.05,]
write.table(sig, file=paste0("./sig_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#