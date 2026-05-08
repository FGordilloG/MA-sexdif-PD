# Single cell RNA-seq data analysis
# Feature Selection
# Autor: Fernando Gordillo Gonzalez
# Args: nombre_estudio

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
dir.create("results/3.Dimentional_Reduction")
dir.create("results/3.Dimentional_Reduction/Features_Selection")
dir.create("figures/3.Dimentional_Reduction")
dir.create("figures/3.Dimentional_Reduction/Features_Selection")

# Load data
sce = readRDS(paste("results/2.Normalization/sce_",tolower(estudio),"_norm.rds", sep=''))

#---------------------------- Feature selection -------------------------------#

# Compute gene variance
dec.sce = modelGeneVar(sce, block=sce$Sample ,BPPARAM=MulticoreParam(20))

# Visualizing the fit:
blocked <- dec.sce$per.block
dir.create("./figures/3.Dimentional_Reduction/Features_Selection/Sample_varianza")
for (i in colnames(blocked)) {
  jpeg(file=paste0("./figures/3.Dimentional_Reduction/Features_Selection/Sample_varianza/var_sample_",i,".jpeg"), width=6, height=4, units="in", res=300)
  current <- blocked[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Media de la expresión log-normalizada", ylab="Varianza log-norm")
  curfit <- metadata(current)
  curve(curfit$trend(x), col='yellow', add=TRUE, lwd=2) 
  dev.off()
}

# Ordering by most biological-variant genes for inspection.
head(arrange(as.data.frame(dec.sce[,-7]),desc(bio)))

# Select the 25% most variable genes for downstream analysis.
chosen = getTopHVGs(dec.sce, prop = 0.25)
write(chosen, file="results/3.Dimentional_Reduction/Features_Selection/HVGs.txt", ncolumns=1)

# Subsetting the HVGs 
sce.hvg = sce[chosen,]
dim(sce.hvg)


# Keep the original SCE with HGV subset
rowSubset(sce, "HVGs") <- chosen

# Plot 25 most expressed genes
jpeg(file="./figures/3.Dimentional_Reduction/Features_Selection/most_expressed_selected.jpeg", width=6, height=4, units="in", res=300)
plotHighestExprs(sce.hvg, n = 25, exprs_values = "counts")+
  theme(legend.position = "none") +
  ggtitle("25 genes más expresados")+
  xlab("Conteos")+ ylab("Gen")
dev.off()

# --------------------------------- Save data --------------------------------#

saveRDS(sce.hvg, paste("results/3.Dimentional_Reduction/Features_Selection/sce_",tolower(estudio),"_feat.select.rds", sep=''))
saveRDS(sce, paste("results/3.Dimentional_Reduction/Features_Selection/sce_",tolower(estudio),"_feat.select_original.rds", sep=''))


# -----------------------------------------------------------------------------#
