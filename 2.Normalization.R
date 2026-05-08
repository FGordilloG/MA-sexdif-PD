# Single cell RNA-seq data analysis
# Normalization
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID Ensembl_ID_plot

#----------------------------- Load packages ----------------------------------#

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(dplyr)
library(Matrix)
library(viridis)
library(BiocParallel)

#---------------------- Make directories & Load data --------------------------#

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
setwd(estudio)
getwd()

# Create directories
dir.create("results/2.Normalization")
dir.create("figures/2.Normalization")

# Load data
sce = readRDS(paste("results/1.QC/sce_",tolower(estudio),"_processed.rds", sep=''))

#---------------------- Normalization by deconvolution ------------------------#

set.seed(789)
# Make cell-clusters for normalization
quick_clusters = quickCluster(sce,BPPARAM=MulticoreParam(6))
# Calculate size factors
sce = computeSumFactors(sce,clusters = quick_clusters,BPPARAM=MulticoreParam(6))
# Scaling and log-transforming
sce = logNormCounts(sce, log=TRUE,BPPARAM=MulticoreParam(6))

print("Resumen de los resultados de normalizacion por deconvolucion:")
summary(sizeFactors(sce))

#---------------------------- Comparative plots -------------------------------#

sce@colData$sum_norm = Matrix::colSums(logcounts(sce))
metadata_norm = as.data.frame(sce@colData)

png(file="./figures/2.Normalization/Size_factors.png", width=6, height=4, units="in", res=300)
hist(sizeFactors(sce), main = "Histograma de los factores de escalado",
     xlab = "Factores de escalado", ylab ="Frecuencia") 
dev.off()

if (!is.na(args[2])){
  jpeg(file="./figures/2.Normalization/Expression_compare_gene.jpeg", width=16, height=6, units="in", res=500)
  gridExtra::grid.arrange(
    plotExpression(sce, x="Sample", exprs_values="counts",
                   colour_by="cond_sex", features=args[2]) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      theme(legend.position = c(0.05,0.85)) +  xlab("Paciente")+
      ggtitle("No normalizado") + scale_color_viridis_d("Condicion"),
    plotExpression(sce, x="Sample", exprs_values="logcounts",
                   colour_by="cond_sex", features=args[2]) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      xlab("Paciente") + theme(legend.position = "none") +
      ggtitle("Normalizado") + scale_color_viridis_d(),
    ncol=2 )
  dev.off()
}
jpeg(file="./figures/2.Normalization/UMIS_comparison_box.jpeg", width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(
metadata_norm %>%
  ggplot(aes(x=Sample, y=sum, fill=cond_sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Conteos por paciente no normalizado")+
  scale_fill_viridis_d(option="turbo"),
metadata_norm %>%
  ggplot(aes(x=Sample, y=sum_norm, fill=cond_sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Conteos por paciente normalizado")+
  scale_fill_viridis_d(option="turbo"),
ncol=2)
dev.off()

jpeg(file="./figures/2.Normalization/UMIS_comparison.jpeg", width=8, height=12, units="in", res=300)
gridExtra::grid.arrange(
plotColData(sce,y = "sum",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Conteos por paciente no normalizado")+
  xlab("Paciente")+
  ylab("Conteos por célula")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d()+
  theme(legend.position = "none"),
plotColData(sce,y = "sum_norm",x = "Sample",colour_by = "cond_sex")+
  ggtitle("Conteos por paciente normalizado")+
  xlab("Paciente")+
  ylab("Conteos por célula")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_color_viridis_d("Condicion") +
  theme(legend.position = "bottom"),
  nrow=2)
dev.off()

# --------------------------------- Save data --------------------------------#

saveRDS(sce, paste("results/2.Normalization/sce_",tolower(estudio),"_norm.rds", sep=''))

# -----------------------------------------------------------------------------#
