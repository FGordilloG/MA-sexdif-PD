# Function for removing batch effect
# Author: Fernando Gordillo Gonzalez

#------------------------------Loading packages--------------------------------#

library(SingleCellExperiment)
library(ggplot2)
library(harmony)
library(scran)
library(scater)

batch_correction = function(sce, vars_corr){
  pcs=reducedDim(sce,"PCA")
  metadata=colData(sce)
  
  harmonized_pcs <- HarmonyMatrix(
    data_mat  = pcs,       # Matrix with coordinates for each cell (row) along many PCs (columns)
    meta_data = metadata,  # Dataframe with information for each cell (row)
    vars_use  = vars_corr, # Column in meta_data that defines dataset for each cell
    do_pca    = FALSE)     # Since we are providing PCs, do not run PCA
     
  
  # attributes(harmonized_pcs) <- list(class = c("reduced.dim.matrix","matrix"), se pierden los atributos pero esto no funciona (borra col/rownames)
  #                                    percentVar = attr(pcs, "percentVar"),
  #                                    varExplained = attr(pcs, "varExplained"),
  #                                    rotation = attr(pcs, "rotation"))
  harmonized_pcs= harmonized_pcs[colnames(sce),]
  reducedDim(sce,"PCA")=harmonized_pcs
  
  return(sce)
}
