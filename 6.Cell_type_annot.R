# Single cell RNA-seq data analysis
# Cell type annotation
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID sce_ref_path genelist_path recap

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
library(SCINA)
library(scmap)
library(SingleR)
library(BRETIGEA)
library(CHETAH)
library(gridExtra)

#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
setwd(estudio)
getwd()

# Create directories
dir.create("results/5.Cell_annotation")
dir.create("figures/5.Cell_annotation")
dir.create("figures/5.Cell_annotation/ind_results")

# Load data
sce = readRDS(paste("results/4.Clustering/sce_",tolower(estudio),"_Clust.rds", sep=''))
rownames(sce)=rowData(sce)$SYMBOL
sce.matrix = as.matrix(logcounts(sce))
# sce.matrix = logcounts(sce)
rownames(sce.matrix)=rowData(sce)$SYMBOL

# Load reference data
sce.ref=readRDS(as.character(args[2]))
geneSets = readRDS(as.character(args[3]))

# Recap annotation
skip=NULL
all_res=data.frame(row.names = colnames(sce))
if (!is.na(args[4])){
  all_res=read.delim("results/5.Cell_annotation/matrix_annot.tsv")
  skip=colnames(all_res)
}

# Plot individual result function
plot_ind <- function(sce, annot,pkg){
  sce$cell_type=annot
  
  jpeg(file=paste0("figures/5.Cell_annotation/ind_results/",pkg,".jpeg"), width=20, height=6, units="in", res=300)
  grid.arrange(
    plotReducedDim(sce, dimred="PCA",colour_by="cell_type", text_by="cluster") +
      ggtitle("Gráfico de PCA por tipo celular") + theme(legend.position = "right") +
      scale_color_viridis_d("Tipo celular",option="H"),
    plotReducedDim(sce, dimred="TSNE_90",colour_by="cell_type", text_by="cluster")+
      ggtitle("TSNE (p=90)") + theme(legend.position = "none") +
      scale_color_viridis_d(option="H"),
    plotReducedDim(sce, dimred="UMAP_300",colour_by="cell_type", text_by="cluster")+
      ggtitle("UMAP (n=300)") + theme(legend.position = "none") +
      scale_color_viridis_d(option="H"),
    nrow=1)
  dev.off()
  
  tab <- table(Assigned=annot, Cluster=sce$cluster)
  jpeg(file=paste0("figures/5.Cell_annotation/ind_results/",pkg,"_heatmap.jpeg"))
    pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
  dev.off()
}

#-------------------------------- Annotation ----------------------------------#
#-------------- SingleR --------------####
# compare the single-cell expression profiles with previously annotated reference datasets
if (is.na(args[4]) | !("SingleR_res" %in% skip)){
cat("\n*SingleR*\n")

SingleR_res = SingleR(test=sce, ref=sce.ref, labels=sce.ref$cell_type_alias_label)

# heatmap:Ideally, each cell should exhibit a high score in one label relative to all
# of the others. Pero en células parecidas o relacionadas es raro.
jpeg(file="figures/5.Cell_annotation/ind_results/SingleR_score.jpeg")
plotScoreHeatmap(SingleR_res)
dev.off()

## Extract and save celltypes
all_res$SingleR_res=SingleR_res$labels
write.table(all_res,file="results/5.Cell_annotation/matrix_annot.tsv", sep="\t", row.names = F)
rm(SingleR_res)

## PLOTS
plot_ind(sce, all_res$SingleR_res, "SingleR")

}

#------------- CHETAH ----------####
if (is.na(args[4]) | !("chetah_res" %in% skip)){
cat("\n*CHETAH*\n")

# SUPPLEMENTARY INFO:
  # Acording to CHETAH developers, it works better with around 300 cells per 
  # cell type, so it is recommended to adapt our reference with something like:
# sce.ref=sce.ref[,c(sample(which(sce.ref$cell_type_alias_label=="Neuron"), 300), which(sce.ref$cell_type_alias_label!="Neuron"))]

# CHETAH
chetah_res <- CHETAHclassifier(input = sce, ref_cells = sce.ref, ref_ct = "cell_type_alias_label",input_c = "logcounts",ref_c = "logcounts")

## Plot the classification
jpeg(file="figures/5.Cell_annotation/ind_results/CHETAH_nodes.jpeg")
PlotCHETAH(chetah_res,redD="TSNE_90", interm=T)
dev.off()

# Filtering results
chetah_res=chetah_res$celltype_CHETAH
chetah_res[chetah_res=="Unassigned"]="unassigned"
chetah_res[grep("Node", chetah_res)]="unassigned"

## Extract and save celltypes
all_res$chetah_res=chetah_res[colnames(sce)]
write.table(all_res,file="results/5.Cell_annotation/matrix_annot.tsv", sep="\t", row.names = F)
rm(chetah_res)

## PLOTS
plot_ind(sce, all_res$chetah_res, "CHETAH")

}
#--------------- SCINA -------------####
if (is.na(args[4]) | !("SCINA_res" %in% skip)){
cat("\n*SCINA*\n")
scina_out=SCINA(sce.matrix, geneSets)
SCINA_res=scina_out$cell_labels
SCINA_res[SCINA_res=="unknown"]="unassigned"

## Extract and save celltypes
all_res$SCINA_res=SCINA_res
write.table(all_res,file="results/5.Cell_annotation/matrix_annot.tsv", sep="\t", row.names = F)
rm(SCINA_res)

## PLOTS
plot_ind(sce, all_res$SCINA_res, "SCINA")
}
#------------- BRETIGEA ------------####
# cuantos más genes marcadores (nMarkers), más probabilidades  de promediar 
# cambios de expresión específicos del tipo celular que pueda ocurrir en los
# grupos de la muestra
if (is.na(args[4]) | !("BRETIGEA_res" %in% skip)){

cat("\n*BRETIGEA*\n")
df_res = brainCells(sce.matrix,species = "human", nMarker = 100, celltypes = c("ast", "mic", "neu", "oli"))
BRETIGEA_res = colnames(df_res)[apply(df_res,1,which.max)]
BRETIGEA_res = as.character(factor(BRETIGEA_res, levels = c("ast","mic","neu","oli"),labels = c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte")))

## Extract and save celltypes
all_res$BRETIGEA_res=BRETIGEA_res
write.table(all_res,file="results/5.Cell_annotation/matrix_annot.tsv", sep="\t", row.names = F)
rm(BRETIGEA_res)

## PLOTS
plot_ind(sce, all_res$BRETIGEA_res, "BRETIGEA")
}

#-------------------- Integration and cell-type inference ---------------------####

all_res$cell_type=apply(all_res,1, function(x) names(which.max(table(as.character(x)))))
all_res$n_suport=apply(all_res[-length(all_res)],1, function(x) max(table(as.character(x))))
all_res$cell_type[which(all_res$n_suport==1)]="unassigned"
apply(all_res, 2, table)

write.table(all_res %>% tibble::rownames_to_column("Barcode"),file="results/5.Cell_annotation/matrix_annot.tsv", sep="\t", row.names = F)


#----------------------------- Cell type plots -------------------------------####
sce$cell_type=all_res$cell_type
jpeg(file="figures/5.Cell_annotation/PCA_cell_type.jpeg", width=8, height=6, units="in", res=300)
plotReducedDim(sce, dimred="PCA",colour_by="cell_type", text_by="cluster")+
  ggtitle("Gráfico de PCA por tipo celular")+
  scale_color_viridis_d("Tipo celular",option="H")
dev.off()


jpeg(file="figures/5.Cell_annotation/TSNE_cell_type.jpeg", width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE_25",colour_by="cell_type", text_by="cluster")+
    ggtitle("TSNE (p=25)")+
    theme(legend.position = "none") +
    scale_color_viridis_d(option="H"),
  plotReducedDim(sce, dimred="TSNE_90",colour_by="cell_type", text_by="cluster")+
    ggtitle("TSNE (p=90)")+
    scale_color_viridis_d("Tipo celular",option="H"), 
  nrow=1)
dev.off()


jpeg(file="figures/5.Cell_annotation/UMAP_cell_type.jpeg", width=16, height=6, units="in", res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP_30",colour_by="cell_type", text_by="cluster")+
    ggtitle("UMAP (n=30)")+
    theme(legend.position = "none") +
    scale_color_viridis_d(option="H"),
  plotReducedDim(sce, dimred="UMAP_300",colour_by="cell_type", text_by="cluster")+
    ggtitle("UMAP (n=300)")+
    scale_color_viridis_d("Tipo celular",option="H"), 
  nrow=1)
dev.off()


# heatmap comparativo clusters-tipos
tab <- table(Assigned=sce$cell_type, Cluster=sce$cluster)
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
jpeg("figures/5.Cell_annotation/heatmap_type_vs_clust.jpeg", width=8, height=6, units="in", res=300)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

# Cell types proportions per patient
tab <- table(Sample=sce$Sample, cell_type=sce$cell_type)

jpeg(file="./figures/5.Cell_annotation/Cell_prop_patients.jpg", width=14, height=6, units="in", res=250)
ggpubr::ggarrange(
  ggplot(as.data.frame(round(tab/rowSums(tab) * 100, 2)), aes(Sample,Freq, fill = cell_type)) + 
    geom_bar(position = "fill",
             stat = "identity", colour="black", size=0.3) +
    ggtitle("Cell proportions per patient") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(aes(label=Freq),position = position_fill(vjust = 0.5),size=3) +
    scale_fill_viridis_d(option="H", begin = 0.1),
  ggplot(as.data.frame(tab), aes(Sample,Freq, fill = cell_type)) + 
    geom_bar(position = "stack",
             stat = "identity", colour="black", size=0.3) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_text(aes(label=Freq),position = position_stack(vjust = 0.5),size=2) +
    ggtitle("Cell numbers by patient") +
    scale_fill_viridis_d(option="H", begin = 0.1),
  nrow=1, common.legend = TRUE, legend="bottom"
)
dev.off()

#------------------------------ Plots per cell --------------------------------####
dir.create("figures/5.Cell_annotation/per_cell")

# Plot cell type
plot_cell <- function(sce, cell,dim){
  sce=sce[,sce$cell_type==cell]
  plot = plotReducedDim(sce, dimred=dim, text_by="cluster") + ggtitle(cell) 
  return(plot)
}

for (red_dim in c("TSNE_90","UMAP_300")){
  plots=list()
  for (cell_type in unique(sce$cell_type)){
    plots[[cell_type]]= plot_cell(sce,cell_type,red_dim)
  }
  jpeg(file=paste0("figures/5.Cell_annotation/per_cell/", red_dim,".jpeg"))
  nRow <- floor(sqrt(length(plots)))
  do.call("grid.arrange", c(plots, nrow=nRow))
  dev.off()
}

#---------------------------------- Save --------------------------------------####

rownames(sce)=rowData(sce)$ENSEMBL
saveRDS(sce, paste("results/5.Cell_annotation/sce_",tolower(estudio),"_annot.rds", sep=''))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#