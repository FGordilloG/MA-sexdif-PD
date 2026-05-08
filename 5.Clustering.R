# Single cell RNA-seq data analysis
# Clustering
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID centers k

#----------------------------- Load packages ----------------------------------####

library(SingleCellExperiment)
library(scater)
library(scran)
library(bluster)
library(uwot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(PCAtools)
library(viridis)
library(BiocParallel)
library(pheatmap)
library(grid)


#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
setwd(estudio)
getwd()

# # Create directories
dir.create("results/4.Clustering")
dir.create("figures/4.Clustering")
dir.create(paste("figures/4.Clustering/",args[2],"_",args[3], sep=''))

# Load data
sce = readRDS(paste("results/3.Dimentional_Reduction/sce_",tolower(estudio),"_Reduct.Dims.rds", sep=''))

#-------------------------------- Clustering ----------------------------------####

## TwoStep clust ####
# k-means is initially used to obtain representative centroids that are subjected to graph-based clustering.
# centers determines the trade-off between speed and fidelity.
set.seed(789)
g.clust = clusterCells(sce, use.dimred="PCA", full=TRUE,
                                BLUSPARAM=TwoStepParam(
                                  first=KmeansParam(centers=as.integer(args[2])),
                                  second=NNGraphParam(k=as.integer(args[3]))))
sce$cluster = g.clust$clusters
write.table(data.frame(barcode=colnames(sce),cluster=sce$cluster),
            file=paste0("./figures/4.Clustering/",args[2],"_",args[3],"/clusters_",estudio,".tsv"),
            sep="\t", row.names = F, col.names = F)
table(sce$cluster)



#----------------------------- Clustering plots -------------------------------#####
jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/PCA_Clusters.jpeg",sep=''), width=8, height=6, units="in", res=300)
plotReducedDim(sce, dimred="PCA",colour_by="cluster", text_by="cluster")+
  ggtitle("Gráfico de PCA por clusters")+
  scale_color_viridis_d("Cluster", option="H")
dev.off()

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/PCA_4_Clusters.jpeg",sep=''), width=16, height=12, units="in", res=300)
plotReducedDim(sce, dimred="PCA",  ncomponents=4, colour_by="cluster")+
  scale_color_viridis_d("Cluster", option="H") + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Grafico 4 primeros PCA por condicion")
dev.off()

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/TSNE_Clusters.jpeg",sep=''), width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE_25",colour_by="cluster", text_by="cluster")+
    ggtitle("TSNE por clusters (p=25)")+
    theme(legend.position = "none") +
    scale_color_viridis_d(option="H"),
  plotReducedDim(sce, dimred="TSNE_90",colour_by="cluster", text_by="cluster")+
    ggtitle("TSNE por clusters (p=90)")+
    scale_color_viridis_d("Cluster", option="H"), 
  nrow=1)
dev.off()

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/UMAP_Clusters.jpeg",sep=''), width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP_30",colour_by="cluster", text_by="cluster")+
    ggtitle("UMAP por clusters (n=30)")+
    theme(legend.position = "none") +
    scale_color_viridis_d(option="H"),
  plotReducedDim(sce, dimred="UMAP_300",colour_by="cluster", text_by="cluster")+
    ggtitle("UMAP por clusters (n=300)")+
    scale_color_viridis_d(option="H"), 
  nrow=1)
dev.off()

#-------------------------------------Save-------------------------------------# ####

saveRDS(sce, paste("results/4.Clustering/sce_",tolower(estudio),"_Clust.rds", sep=''))

#--------------------------- Clustering evaluation ----------------------------####
# Poorly separated clusters can be prioritized for manual inspection to determine the relevance 
# of the differences to its neighboring clusters, while heterogeneous clusters may be subject to
# further clustering to identify internal structure.

# Silhouette width ##
# For each cell, compute the average distance to all cells in the same cluster and the average
# distance to all cells in another cluster, taking the minimum. it naturally captures both
# underclustering and overclustering, we can obtain a “sensible” initial choice for the number
# of clusters by maximizing the average silhouette width

silhouette = as.data.frame(approxSilhouette(reducedDim(sce, "PCA"), clusters=sce$cluster))
silhouette$closest = factor(ifelse(silhouette$width > 0, sce$cluster, silhouette$other))
silhouette$cluster = sce$cluster

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/Silhouette_width.jpeg", sep=''), width=10, height=6, units="in", res=300)
ggplot(silhouette, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  ggtitle("Ancho de silueta por cluster") +
  xlab("Grupos") +
  ylab("Ancho") +
  labs(color="Cluster más\ncercano")
dev.off()


## Purity ##
# defined for each cell as the proportion of neighboring cells that are assigned to the same cluster.

purity = as.data.frame(neighborPurity(reducedDim(sce, "PCA"), sce$cluster))
purity$cluster = sce$cluster

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/Purity.jpeg", sep=''), width=10, height=6, units="in", res=300)
ggplot(purity, aes(x=cluster, y=purity, colour=maximum)) +
  ggbeeswarm::geom_quasirandom() +
  ggtitle("Pureza por cluster") +
  xlab("Grupos") +
  ylab("Pureza") +
  labs(color="Cluster\nmayoritario")
dev.off()
# clusters with low purity values may warrant more careful inspection (closely related subpopulations)

## Modularity ##
# Separation between communities. Larger modularity values indicate that there most edges occur
# within clusters, suggesting that the clusters are sufficiently well separated to avoid edges
# forming between neighboring cells in different clusters

ratio = pairwiseModularity(g.clust$objects$second$graph, g.clust$objects$second$communities$membership, as.ratio=TRUE)

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/Pheatmap_Modularity.jpeg",sep=''), width=10, height=6, units="in", res=300)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE, col=viridis::turbo(100), angle_col = "0", main="Modularity")
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", y=-0.02, gp=gpar(fontsize=10))
dev.off()


cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1),
                                                  mode="upper", weighted=TRUE, diag=FALSE)
set.seed(789)
jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/Graph_Modularity.jpeg",sep=''), width=8, height=8, units="in", res=300)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5,
            layout=igraph::layout_with_lgl, vertex.label.font=2, vertex.label.cex=1.7,
            vertex.label.color="black",  vertex.color=turbo(8,begin = 0.1))
dev.off()
rm(cluster.gr)


## Cluster stability ##
# stable to perturbations to the input data.greater stability also increases the
# likelihood that those conclusions can be reproduced in an independent replicate study

fun_cluster <- function(x) {
  clusterCells(sce, use.dimred="PCA",
               BLUSPARAM=TwoStepParam(
                 first=KmeansParam(centers=100),
                 second=NNGraphParam(k=5)))
}

set.seed(789)
ratios <- bootstrapStability(reducedDim(sce,"PCA"), FUN=fun_cluster, clusters=sce$cluster)

jpeg(file=paste("./figures/4.Clustering/",args[2],"_",args[3],"/Pheatmap_bootstrap.jpeg",sep=''), width=10, height=8, units="in", res=300)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95,
            height=0.95, name="vp", just=c("right","top"))), action="prepend")
pheatmap(ratios, cluster_cols=FALSE, cluster_rows=FALSE, col=viridis::turbo(100),
         angle_col = "0", main="Cluster stability")
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", y=-0.02, gp=gpar(fontsize=10))
dev.off()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#