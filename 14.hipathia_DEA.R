# sn-RNA-seq sex-difs atlas in PD
# Signaling Pathway analysis - DEA of hipathia (MAST)
# Autor: Fernando Gordillo Gonzalez
# Args: nombre_estudio cell_type

#----------------------------- Load packages ----------------------------------####

library(hipathia)
library(SingleCellExperiment)
library(MAST)
library(edgeR)

#------------------------------- Functions ---------------------------------####

zz_process = function(zlm, sce) {
  
  results = cbind(data.frame(rownames(zlmCond_cells)), data.frame(zlmCond_cells[,3,c(1,3)]))
  colnames(results) = c("feat_ID","lambda","p.value")
  results$feat_name=rowData(sce)[match(results$feat_ID, rownames(sce)), "path.name"]
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
  
  
  control=Hypothesis(control, colnames(zlm@LMlike@modelMatrix))
  case=Hypothesis(case, colnames(zlm@LMlike@modelMatrix))
  
  fc=getLogFC(zlm, contrast0, contrast1)
  sum =summary(zlm, logFC=fc)$datatable
  sum=sum[contrast == "case" & component=='logFC',-c(2,3)]
  sum$SE <- (sum[, "ci.hi"] - sum[, "ci.lo"])/ 3.92
  colnames(sum)[c(1,4)]=c("feat_ID","logFC")
  sum[is.na(sum)]=0
  
  results = merge(results, sum, by ="feat_ID", sort=F)
  
  return(results)
}

#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
cell_type = as.character(args[2])

# Create directories
dir.create(paste0(estudio,"/results/13.hipathia/DEA_MAST/",cell_type), recursive = T)
setwd(paste0(estudio,"/results/13.hipathia/DEA_MAST/",cell_type))
getwd()

results = readRDS(paste0("../../hipathia_",estudio,".rds")) 
dim(assay(results))
assay(results)[1:4, 1:4] # each row an effector subpathway and each column a cell

# SummarizedExperiment of effector pathways of celltype
path_vals <- get_paths_data(results) 
dim(assay(path_vals)) # each row an effector subpathway (1876) and each column a cell
path_vals <- path_vals[,path_vals$cell_type == cell_type]

print(paste0("Number of ", cell_type, ": ", ncol(assay(path_vals)), collapse = " "))

#------------------------------ MAST analysis ---------------------------------####

# 1. Hurdle model
# Construct sca object
sca.cell = SceToSingleCellAssay(as(path_vals, "SingleCellExperiment"))
rowData(sca.cell)$primerid = rownames(path_vals)

# Model
# zlmCond = zlm(~0+cond_sex+Sample, sca.cell, method = "bayesglm", ebayes = TRUE, parallel = TRUE)
# saveRDS(zlmCond, "zlm_model.rds")
zlmCond <- readRDS("zlm_model.rds")

### CaseFemale - ControlFemale ####
# Test hypothesis
# zlmCond_cells <- lrTest(zlmCond, Hypothesis("cond_sexPD_F - cond_sexC_F"))
# saveRDS(zlmCond_cells, "zlm_lrt_F.rds")
zlmCond_cells <- readRDS("zlm_lrt_F.rds")
  
results = zz_process(zlmCond_cells,sca.cell)
write.table(results, file=paste0("./all_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond,results,"cond_sexPD_F", "cond_sexC_F")
results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resF.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)


### CaseMale - ControlMale ####
# Test hypothesis
zlmCond_cells <- lrTest(zlmCond, Hypothesis("cond_sexPD_M - cond_sexC_M"))
saveRDS(zlmCond_cells, "zlm_lrt_M.rds")

results = zz_process(zlmCond_cells,sca.cell)
write.table(results, file=paste0("./all_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond,results,"cond_sexPD_M", "cond_sexC_M")
results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)


### (CaseFemale - ControlFemale) - (CaseMale - ControlMale) ####
# Test hypothesis
zlmCond_cells <- lrTest(zlmCond, Hypothesis("(cond_sexPD_F - cond_sexC_F) - (cond_sexPD_M - cond_sexC_M)"))
saveRDS(zlmCond_cells, "zlm_lrt_FM.rds")

results = zz_process(zlmCond_cells,sca.cell)
write.table(results, file=paste0("./all_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)

results= logFC_mast(zlmCond, results, c("cond_sexC_M","cond_sexPD_F"), c("cond_sexC_F","cond_sexPD_M"))
results = results[order(results$p.value), ]
write.table(results, file=paste0("./all_resFM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t",quote = FALSE)
