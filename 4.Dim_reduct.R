# Single cell RNA-seq data analysis
# Dimention reduction
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID batch_correction covariates(,separed)

#----------------------------- Load packages ----------------------------------#

library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(PCAtools)
library(viridis)
library(BiocParallel)

#---------------------- Make directories & Load data --------------------------#

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
setwd(estudio)
getwd()

# Load data
sce = readRDS(paste("results/3.Dimentional_Reduction/Features_Selection/sce_",tolower(estudio),"_feat.select_original.rds", sep=''))


#------------------------------------ PCA ----------------------------------####
set.seed(789)
sce=fixedPCA(sce, subset.row=rowSubset(sce,'HVGs'), rank = 30)

# Keep only PCs before elbow ####
pve = attr(reducedDim(sce, "PCA"), "percentVar")
cum = cumsum(pve)
elbow = PCAtools::findElbowPoint(pve)
print(paste("Número de PCAs elegidos por el corte del codo:",elbow))

## elbow plot
jpeg(file="./figures/3.Dimentional_Reduction/PCA_Variance.jpeg", width=14, height=6, units="in", res=300)
par(mfrow = c(1, 2))
plot(pve, xlab="PC", ylab="Variance(%)", main= "Variance explained", type = "b")
points(x=elbow, y=pve[elbow], pch=16, cex=1.3, col="red")
plot(cum, xlab = "PC", ylab="Variance(%)", main= "Cumulative variance explained",
     ylim = c(0, 100), type = "b")
points(x=elbow, y=cum[elbow], pch=16, cex=1.3, col="red")
text(elbow, cum[elbow],
     labels = "Elbow", cex = 0.6, pos = 3, col = "black")
text(elbow, (cum[elbow]-3),
     labels = paste("Tot. var. expl : ",format(cum[elbow],2)," %", sep=""), cex = 0.6, pos = 4, col = "black")
dev.off()

reducedDim(sce, "PCA") = reducedDim(sce, "PCA")[,1:elbow]

## PLOTS ####
jpeg(file="./figures/3.Dimentional_Reduction/PCA_Condition.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce, dimred="PCA", colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Gráfico PCA por condicion")
dev.off()

#### Batch correction ####
if (!is.na(args[2])){
  source("./Batch_correction.R")
  covariates=strsplit(args[3], ",")[[1]]
  sce=batch_correction(sce,covariates)
}

jpeg(file="./figures/3.Dimentional_Reduction/PCA_Patient.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce, dimred="PCA", colour_by="Sample")+
  scale_color_viridis_d("Sample") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Gráfico PCA por paciente")
dev.off()

jpeg(file="./figures/3.Dimentional_Reduction/PCA_4_Condition.jpeg", width=16, height=12, units="in", res=300)
plotReducedDim(sce, dimred="PCA",  ncomponents=4, colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Grafico 4 primeros PCA por condicion")
dev.off()

#----------------------------------- t-SNE ------------------------------------#

set.seed(789)
sce = runTSNE(sce, dimred="PCA", perplexity=25, name = "TSNE_25")
tsne25=plotReducedDim(sce, dimred="TSNE_25", colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") + theme(legend.position = "none")+
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("perplexity=25")

set.seed(789)
sce = runTSNE(sce, dimred="PCA", perplexity=90, name = "TSNE_90")
tsne90=plotReducedDim(sce, dimred="TSNE_90", colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") +
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("perplexity=90")

jpeg(file="./figures/3.Dimentional_Reduction/TSNE_Condition.jpeg", width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(tsne25,tsne90,nrow=1)
dev.off()

#----------------------------------- UMAP -------------------------------------#

set.seed(789)
sce = runUMAP(sce, dimred="PCA", n_neighbors= 30, name = "UMAP_30")
umap30= plotReducedDim(sce, dimred="UMAP_30", colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") + theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("n_neighbors= 30")


sce = runUMAP(sce, dimred="PCA",  n_neighbors= 300, name = "UMAP_300")
umap300= plotReducedDim(sce, dimred="UMAP_300", colour_by="cond_sex")+
  scale_color_viridis_d("Condicion") + 
  theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("n_neighbors= 300")


jpeg(file="./figures/3.Dimentional_Reduction/UMAP_Condition.jpeg", width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(umap30,umap300,nrow=1)
dev.off()

# --------------------------------- Save data --------------------------------#

saveRDS(sce, paste("results/3.Dimentional_Reduction/sce_",tolower(estudio),"_Reduct.Dims.rds", sep=''))

# -----------------------------------------------------------------------------#
