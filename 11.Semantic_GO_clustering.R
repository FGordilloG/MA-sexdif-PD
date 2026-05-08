# sn-RNA-seq sex-difs atlas in PD
# Semantic clustering among GO-sets
# Autor: Fernando Gordillo Gonzalez
# Args: cell_type tissue

#----------------------------- Load packages ----------------------------------####

library(stringr)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(dplyr)

#---------------------- Make directories & Load data -----------------------####

setwd("~/PDsex_atlas/Integracion/MA")

## Global variables
# Generic words to avoid in the word-clouds
no_words=c("regulation", "process","cell", "signaling", "pathway","particle",
           "nuclear","negative", "positive","tissue", "activity","cardiac","muscle")
medida = "Rel"
organismo = "org.Hs.eg.db"
algoritmo = "personalizado"

#---------------------- Make directories & Load data -----------------------####
## BETWEEN TISSUES (same cells) --------------####
cell="Oligodendrocyte"
labels=c("up","up") #(1 = SN ; 2 = C)
pvalor = 0.05

weight01=list()
a = read.delim(paste0("SN/",cell,"/Enrichment_weight01_",cell,"_",labels[1],".tsv"))
weight01[["SN"]] = a[a$p.value<pvalor,"GO.ID"]
a = read.delim(paste0("C/",cell,"/Enrichment_weight01_",cell,"_",labels[2],".tsv"))
weight01[["C"]] = a[a$p.value<pvalor,"GO.ID"]

# Unique terms
all.terms = unique(unname(unlist(weight01)))

# heatmap matrix
m = matrix(0, nrow = length(all.terms), ncol = 2)
rownames(m) = all.terms
colnames(m) = c("SN","C")

for (tiss in colnames(m)) {
  m[rownames(m) %in% weight01[[tiss]],tiss]=1
}


## BETWEEN CELLS (same tissue and sex) --------------####
# t="SN"
# sex="male"
# label=ifelse(sex=="female","up","down") #(UP = FEMALE ; DOWN = MALE)
# # cell_types=c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte")
# cell_types=c("Astrocyte", "Microglia")
# pvalor = 0.05
# 
# weight01=list()
# for (cell in cell_types) {
#   a = read.delim(paste0(t,"/",cell,"/Enrichment_weight01_",cell,"_",label,".tsv"))
#   weight01[[cell]] = a[a$p.value<pvalor,"GO.ID"]
# }
# # Unique terms
# all.terms = unique(unname(unlist(weight01)))
# 
# # heatmap matrix
# m = matrix(0, nrow = length(all.terms), ncol = length(cell_types))
# rownames(m) = all.terms
# colnames(m) = cell_types
# 
# for (cell in cell_types) {
#   m[rownames(m) %in% weight01[[cell]],cell]=1
# }

## ---------------------- GO semantic clustering ---------------------------####
# 4. Similarity matrix

mat = GO_similarity(all.terms, ont = "BP", db = organismo , measure = medida)

# 4. Clustering
heatp = Heatmap(m, col = c("white", "#632C63"), cluster_columns = FALSE, cluster_rows = FALSE, 
                column_names_gp = gpar(fontsize = 9),
                show_heatmap_legend = T,
                heatmap_width = unit(1, "cm"),
                heatmap_legend_param = list(
                                      title = "GO ID",
                                      labels = c("Significant", " Not significant"))
                )

register_clustering_methods(
  personalizado = function(mat) binary_cut(mat, try_all_partition_fun = TRUE, cutoff = 0.85)
)
set.seed(789)
df_fisher = simplifyGO(mat, method = "personalizado" ,plot=F)
# df_fisher = simplifyGO(mat, method = "louvain" ,plot=F)
df_fisher = cbind(df_fisher, m[df_fisher$id,])

# PLOT TISSUES ####
# Para anotar solo los clusters que tengan mezcla
go_list=df_fisher %>% group_by(cluster) %>% summarise(go_terms = list(id)) 
go_list=structure(pull(go_list,go_terms),names=pull(go_list,cluster))

cluster_word = as.integer(df_fisher %>% group_by(cluster) %>% summarise(SN = sum(SN),C = sum(C)) %>% filter(SN>1 & C>1) %>% pull(cluster))

align_to=as.factor(df_fisher$cluster)
align_to[!align_to %in% cluster_word]=NA
align_to=droplevels(align_to)

ht=ht_clusters(mat, df_fisher$cluster, col = c("white", "#632C63"), draw_word_cloud = F, 
            order_by_size = TRUE, min_term = 3, column_title = paste("GO terms over-represented in",cell,paste(labels, collapse = " ")),
            use_raster = FALSE, ht_list = heatp, run_draw = F)

ht=ht[,1:2] + 
  rowAnnotation(go = anno_word_cloud_from_GO(align_to, go_list, 
                                             fontsize_range = c(7, 10), max_words=6,
                                             word_cloud_grob_param=list(word_space = unit(10, "mm"),
                                                                        line_space = unit(8, "pt"),
                                                                        max_width = unit(70, "mm")),
                                             bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA"),
                                             exclude_words=c(no_words,"reninangiotensin","presure", "blood","arterial"))
                )

png(paste0("./htclust_tissue/",cell,"_SN-",labels[1],"_C-",labels[2],"_comunes_SN4.png"), width = 9, height = 6, units="in",res = 300)
draw(ht, main_heatmap = "Similarity", gap=unit.c(unit(4,"pt"),unit(1,"pt")))
dev.off()


# PLOT CELLS ####
# dir.create(paste(t,"htclust_combi",sep="/"))
setwd(paste(t,"htclust_combi",sep="/"))

if(length(cell_types)==4){
  png(paste0("./all_cells_",sex,"_SN.png"), width = 2780, height = 1480, res = 200)
} else {
  png(paste0("./",paste(substr(cell_types,1,3), collapse = "_"),"_",sex,"_SN.png"), width = 3580, height = 2080, res = 300)
}
ht_clusters(mat, df_fisher$cluster, col = c("white", "#632C63"), draw_word_cloud = TRUE,  fontsize_range = c(8, 11), max_words = 6,
            order_by_size = TRUE, bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA"),
            min_term = 4,
            column_title = paste("GO terms over-represented in",sex,"PD"),
            use_raster = FALSE, ht_list = heatp, exclude_words=c(no_words))
dev.off()

