# Single cell RNA-seq data analysis
# Quality Control (QC)
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID ncells UMIs ngenes %mito

#----------------------------- Load packages ----------------------------------#

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)
library(viridis)
library(BiocParallel)
library(scDblFinder)

#---------------------- Make directories & Load data --------------------------#

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
setwd(estudio)
getwd()

# Create directories
dir.create("results/1.QC")
dir.create("figures")
dir.create("figures/1.QC")

# Load data
sce = readRDS(paste("results/0.Preprocesamiento/sce_",tolower(estudio),"_raw.rds", sep=''))

#------------------------------------ QC -----------------------------------####
cat("\nDimensiones originales\n")
dim(sce)

# Creating mitocondrial genes subset
mito_genes=rowData(sce)[grep("^MT-",rowData(sce)$SYMBOL),"ENSEMBL"]
cat("\nMitochondrial genes:\n")
length(mito_genes)

# QC metrics
sce <- addPerCellQC(sce, subsets=list(Mito= mito_genes),BPPARAM=MulticoreParam(10))

#--------------------------- Filtering features ----------------------------####

# nCells for each gene
rowData(sce)$nCells <- Matrix::rowSums(counts(sce)>0)

## Non-coding genes (https://www.genenames.org/about/guidelines/) ####
# tRNA
trna= rownames(sce)[grep("^TR[A-Z]-",rowData(sce)$SYMBOL)]
print("Number of tRNA:")
length(trna)

# miRNA
mirna=rownames(sce)[grep("^(MIR|mir-)",rowData(sce)$SYMBOL)]
print("Number of miRNA:")
length(mirna)

# lncRNA
lncrna=rownames(sce)[grep("(^LINC[0-9]{5}|-(AS|[DIO]T)[0-9]?$)",rowData(sce)$SYMBOL)]
print("Number of lncRNA:")
length(lncrna)

# Small nuclear RNAs
smrna=rownames(sce)[grep("^RNU",rowData(sce)$SYMBOL)]
print("Number of smRNA:")
length(smrna)

# Small nucleolar RNAs
snorna=rownames(sce)[grep("^(SNOR[DA]|SCARNA)",rowData(sce)$SYMBOL)]
print("Number of snoRNA:")
length(snorna)

# Ribosomal RNAs
rrna=rownames(sce)[grep("^RNA",rowData(sce)$SYMBOL)]
print("Number of rRNA:")
length(rrna)

# Pseudogenes https://www.biostars.org/p/51456/
psgen = rownames(sce)[grep("^RP[0-9]{2}-",rowData(sce)$SYMBOL)]
print("Number of pseudogenes:")
length(psgen)
#####

low_features = rowData(sce)$nCells > as.integer(args[2])
print(paste("nCells per gene threshold: ",as.character(args[2])))

nc_genes <- c(trna, mirna, lncrna, smrna, snorna, rrna, mito_genes, psgen)
pseudog = rownames(sce) %in% nc_genes
nkeep = !pseudog & low_features
sce = sce[nkeep,]
print("Genes filtrados")
sum(!nkeep)

# Doublet identification
set.seed(789)
sce <- scDblFinder(sce,samples=sce$Sample, BPPARAM=MulticoreParam(10))
print("Detected doublets:")
table(sce@colData$scDblFinder.class)

#---------------------------- Plots before filtering -----------------------####

metadata_bf = as.data.frame(sce@colData)
dir.create("./figures/1.QC/Before_filtering")

cat("\nResumen UMIs detectados:\n")
summary(sce$sum)

cat("\nResumen genes detectados:\n")
summary(sce$detected)

# para saber la media de genes detectados por célula
jpeg(file="./figures/1.QC/Before_filtering/density_plot_features.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce$detected), xlab="Genes detectados", ylab="Densidad", main="Genes detectados por célula")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/NCells_sample.jpeg", width=8, height=6, units="in", res=300)
metadata_bf %>%
  ggplot(aes(x=Sample, fill=cond_sex)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  guides(fill=guide_legend(title="Condicion")) +
  ggtitle("Número de células por paciente")+
  xlab("Paciente")+
  ylab("Número de céĺulas")+
  scale_color_viridis_d(option="viridis")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/Features_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "detected",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Genes detectados por paciente")+
  xlab("Paciente")+
  ylab("Genes por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/Count_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "sum",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Número de conteos por paciente")+
  xlab("Paciente")+
  ylab("Conteos por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/mito_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "subsets_Mito_percent",x = "Sample",colour_by = "cond_sex")+
  ggtitle("% mitocondriales por paciente")+
  xlab("Paciente")+
  ylab("%mito por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/sum_vs_detected.jpeg", width=16, height=6, units="in", res=300)
plotColData(sce, x = "sum", y="detected", colour_by="Sample", other_fields ="condition") +
  facet_wrap(~condition) +
  ggtitle("Counts vs Detected genes")+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  xlab("Counts")+
  ylab("Detected genes")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/sum_vs_mito.jpeg", width=16, height=6, units="in", res=300)
plotColData(sce, x = "sum", y="subsets_Mito_percent", colour_by="Sample", other_fields ="condition") +
  facet_wrap(~condition) +
  ggtitle("Counts vs % mito")+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  xlab("Counts")+
  ylab("% mitochondrial genes")
dev.off()

jpeg(file="./figures/1.QC/Before_filtering/doublets_patient.jpeg", width=8, height=6, units="in", res=300)
ggplot(as.data.frame(table(sce$scDblFinder.class, sce$Sample)), aes(x=Var2, y = Freq, fill=Var1)) +
  geom_bar(stat="identity") +
  scale_fill_manual("Clasificacion",values = c("lightblue","purple"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Identificación de dobletes")+
  xlab("Paciente")+
  ylab("Número de dobletes")
dev.off()

#---------------------- Applying superior limit (3 MADs) -------------------####
# All
sizelib.sup <- isOutlier(sce$sum, log=TRUE, nmads = 3, type="higher")
print("UMIs superior threshold (3 MADS):")
print(attr(sizelib.sup, "thresholds"))
# print(paste("Discarted with 3 MADS:", sum(isOutlier(sce$sum, log=TRUE, nmads = 3, type="higher"))))
#----------------------- Filtering with fixed thresholds ----------------------#
  cat("\n\nFixed thresholds:\n")
  print(paste("UMIs (sum) threshold: ",as.character(args[3])))
  sizelib.ft <- sce@colData$sum < as.integer(args[3])
  print(paste("Feature detected threshold: ",as.character(args[4])))
  nexprs.ft <- sce@colData$detected < as.integer(args[4])
  print(paste("% mito genes threshold: ",as.character(args[5])))
  rmito.ft <- sce@colData$subsets_Mito_percent > as.integer(args[5])

  doublet <- sce@colData$scDblFinder.class == "doublet"
  
  discard.ft <- sizelib.ft | nexprs.ft | rmito.ft  | doublet | sizelib.sup
  
  cat("\nSummary fixed thresholds\n\n")
  print(data.frame( Libsize=sum(sizelib.ft), NExprs=sum(nexprs.ft), MitoProp=sum(rmito.ft), Doublet = sum(doublet), Libsize.sup=sum(sizelib.sup), Total=sum(discard.ft)))
  
  sce_gfilt=sce
  colData(sce_gfilt)$discard = discard.ft
  sce=sce[,!discard.ft]

#----------------------------- Plots after filtering --------------------------#

dir.create("./figures/1.QC/After_filtering")

# Discard plots (with discarded cells coloured)

jpeg(file="./figures/1.QC/After_filtering/Discard_plots.jpeg", width=8, height=24, units="in", res=500)
gridExtra::grid.arrange(
  plotColData(sce_gfilt, x="Sample", y="sum", colour_by="discard", other_fields="condition") + 
    facet_grid(cols=vars(condition), scales = "free_x") +
    scale_color_manual(values=c("grey69","red"))+
    theme(legend.position = "none")+
    theme(plot.title = element_text(hjust=0.5)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Conteos por paciente") +xlab("Paciente") +ylab("Conteos por célula") ,
  plotColData(sce_gfilt, x="Sample", y="detected", colour_by="discard", other_fields="condition") + 
    facet_grid(cols=vars(condition), scales = "free_x") +
    scale_color_manual(values=c("grey69","red"))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) + 
    ggtitle("Genes por paciente") +xlab("Paciente") +ylab("Genes por célula") ,
  plotColData(sce_gfilt, x="Sample", y="subsets_Mito_percent", colour_by="discard", other_fields="condition") + 
    facet_grid(cols=vars(condition), scales = "free_x") +
    scale_color_manual(values=c("grey69","red"))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) + 
    ggtitle("%mito por paciente") +xlab("Paciente") +ylab("%mito por célula"),
  ggplot(as.data.frame(sce_gfilt@colData)[which(sce_gfilt@colData$discard),],
         aes(x=Sample, fill=cond_sex)) + geom_bar() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(fill=guide_legend(title="Condicion")) +
    ggtitle("Células descartadas por paciente")+
    xlab("Paciente")+
    ylab("N descartes")+
    scale_color_viridis_d(option="D")
, nrow = 4)
dev.off()


## Discard cells ##
cat("\n\n\n\nCélulas eliminadas:\n")
sum(sce_gfilt@colData$discard)
cat("\n\n")
rm(sce_gfilt)

cat("\n\n\nDimensiones despues del filtrado total\n")
dim(sce)
cat("\n\n\n")

# Final plots (without discarded)
metadata_af = as.data.frame(sce@colData)

cat("\nResumen UMIs detectados:\n")
summary(sce$sum)

cat("\nResumen genes detectados:\n")
summary(sce$detected)

jpeg(file="./figures/1.QC/After_filtering/density_plot_features.jpeg", width=6, height=4, units="in", res=300)
plot(density(sce$detected), xlab="Genes detectados", ylab="Densidad", main="Genes detectados por célula")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/NCells_sample.jpeg", width=8, height=6, units="in", res=300)
metadata_af %>%
  ggplot(aes(x=Sample, fill=cond_sex)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  guides(fill=guide_legend(title="Condicion")) +
  ggtitle("Número de células por paciente")+
  xlab("Paciente")+
  ylab("Número de céĺulas")+
  scale_color_viridis_d(option="viridis")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/Features_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "detected",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Genes detectados por paciente")+
  xlab("Paciente")+
  ylab("Genes por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/Count_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "sum",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Número de conteos por paciente")+
  xlab("Paciente")+
  ylab("Conteos por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/mito_patient.jpeg", width=8, height=6, units="in", res=300)
plotColData(sce,y = "subsets_Mito_percent",x = "Sample",colour_by = "cond_sex")+
  ggtitle("% mitocondriales por paciente")+
  xlab("Paciente")+
  ylab("%mito por célula")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion")
dev.off()

jpeg(file="./figures/1.QC/After_filtering/sum_vs_detected.jpeg", width=16, height=6, units="in", res=300)
plotColData(sce, x = "sum", y="detected", colour_by="Sample", other_fields ="condition") +
  ggtitle("Counts vs Detected genes")+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  xlab("Counts")+
  ylab("Detected genes") +
  facet_wrap(~condition)
dev.off()

jpeg(file="./figures/1.QC/After_filtering/sum_vs_detected_log10.jpeg", width=16, height=6, units="in", res=300)
plotColData(sce, x = "sum", y="detected", colour_by="Sample", other_fields ="condition") +
  ggtitle("Counts vs Detected genes")+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  xlab("Counts")+
  ylab("Detected genes") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~condition)
dev.off()


jpeg(file="./figures/1.QC/After_filtering/sum_vs_mito.jpeg", width=16, height=6, units="in", res=300)
plotColData(sce, x = "sum", y="subsets_Mito_percent", colour_by="Sample", other_fields ="condition") +
  facet_wrap(~condition) +
  ggtitle("Counts vs % mito")+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  xlab("Counts")+
  ylab("% mitochondrial genes")
dev.off()

# Complexity per cell (https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html)
metadata_af$log10GenesPerUMI <- log10(metadata_af$detected) / log10(metadata_af$sum)
jpeg(file="./figures/1.QC/After_filtering/Complexity_per_cell.jpeg", width=8, height=6, units="in", res=300)
metadata_af %>% ggplot(aes(x=log10GenesPerUMI, color = Sample, fill=Sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

# ---------------- Relative expression of each gene per cell ------------------#

jpeg(file="./figures/1.QC/most_expressed.jpeg", width=6, height=4, units="in", res=300)
plotHighestExprs(sce)
dev.off()

print(paste("Realizado QC del estudio",estudio))

# --------------------------------- Save data --------------------------------#

saveRDS(sce, paste("results/1.QC/sce_",tolower(estudio),"_processed.rds", sep=''))

# -----------------------------------------------------------------------------#
