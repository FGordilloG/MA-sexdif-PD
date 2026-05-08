# sn-RNA-seq sex-difs atlas in PD
# Signaling Pathway analysis - get activation values
# Autor: Fernando Gordillo Gonzalez
# Args: study_ID

#------------------------------ Load packages ------------------------------####
library(hipathia)
library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(Matrix)
library(clusterProfiler)
library("org.Hs.eg.db")

#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])

# Create directories
dir.create(paste0(estudio,"/results/13.hipathia/"), recursive = T)
setwd(paste0(estudio,"/results/13.hipathia/"))
getwd()

# Load data
sce.n = readRDS(paste0("../5.Cell_annotation/sce_",tolower(estudio),"_annot_filt.rds"))

# Pathways
pathways <- load_pathways(species = "hsa")

#---------------------------- Hipathia analysis ------------------------------####
# hipathia assumes that data has been already normalized for correcting any possible 
# sequencing bias (which includes also batch effect correction) 

# 3.1. Translate data to ENTREZID
sce = as(sce.n, "SummarizedExperiment")
rownames(sce) = rownames(sce.n)
rowData(sce) = rowData(sce.n)
assay(sce, "logcounts") = as.matrix(assay(sce.n, "logcounts"))
rm(sce.n)

eg = bitr(rownames(sce), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
sce=sce[eg$ENSEMBL,]
rowData(sce)$ENTREZ = eg$ENTREZID
rownames(sce)=rowData(sce)$ENTREZ

# 3.2. Scaled data
sce_scaled = normalize_data(sce, sel_assay = "logcounts")

# 3.3. Hipathia analysis
results <- hipathia(sce_scaled, pathways, decompose = FALSE, verbose=TRUE)
saveRDS(results, file=paste0("./hipathia_",estudio,".rds"))

