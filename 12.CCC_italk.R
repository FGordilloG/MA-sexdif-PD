# sn-RNA-seq sex-difs atlas in PD
# Cell-Cell Communication analysis
# Autor: Fernando Gordillo Gonzalez
# Args: cell_type contrast


#----------------------------- Load packages ----------------------------------####

library(iTALK)
library(dplyr)
library(tidyr)
library(OmnipathR)
library(ggpubr)
library(viridis)
library(purrr)

#---------------------- Make directories & Load data -----------------------####
# INTERCELL NETWORK
icn <- import_intercell_network(entity_types = "protein")
db_italk= icn %>%
  filter(!(grepl("_",source_genesymbol) | grepl("_",target_genesymbol))) %>%
  mutate(Pair.Name=paste(source_genesymbol,target_genesymbol, sep="_"),
         Classification="intercell") %>%
  distinct(Pair.Name,.keep_all = T) %>%
  rename(Ligand.ApprovedSymbol=source_genesymbol,
         Receptor.ApprovedSymbol=target_genesymbol)


setwd("~/PDsex_atlas/Integracion/MA")
getwd()
dir.create(paste0(tissue, "/Comunicacion_celular"))

tissue= "SN"
sex="male"
umbral=0
ts= "Substantia nigra"
cell_types=c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte")

#---------------------- Function for annotation -----------------------####

anot_path <- function(res,annotation){
  res_annot <- annotated_network(tibble(res), annotation) 
  if (any(!is.na(res_annot$pathway_source))){
    ligand_pth <- aggregate(pathway_source ~ ligand, unique(res_annot[,c("pathway_source","ligand")]), paste, collapse = ", ")
    res[,paste0(annotation,"_ligand")]=ligand_pth[match(res$ligand, ligand_pth$ligand),"pathway_source"]
    res_annot$pathway_source = NULL
  } else { res_annot$pathway_source = NULL }
  
  if (any(!is.na(res_annot$pathway_target))){
    receptor_pth <- aggregate(pathway_target ~ receptor, unique(res_annot[,c("pathway_target","receptor")]), paste, collapse = ", ")
    res[,paste0(annotation,"_receptor")]=receptor_pth[match(res$receptor, receptor_pth$receptor),"pathway_target"]
    res_annot$pathway_source = NULL
  } else { res_annot$pathway_target = NULL }
  return(res)
}

#---------------------- CCC -----------------------####
cell_col<-structure(viridis::turbo(length(cell_types), begin = 0.1),names=cell_types)
all_degs=c()
cell_res=list()
for (cell in cell_types) {
  res=read.delim(paste("../",tissue,cell,"FM/all.genes.tsv",sep="/"))[,c(2,4,6,13)]
  if (sex == "male"){
    res= res %>% filter((logFC < -(umbral) & p.adjust.fdr < 0.05) | (logFC < -(1.5) & p.adjust.fdr < 0.1))
  } else{
    res= res %>% filter((logFC > (umbral) & p.adjust.fdr < 0.05) | (logFC > 1.5 & p.adjust.fdr < 0.1))
  }
  res = res %>% mutate(cell_type=cell, .before=logFC)
  colnames(res)=c("gene","cell_type","logFC","p.value","q.value")
  cell_res[[cell]]=res
  
  all_degs=c(all_degs,res$gene)
}

from=NULL
to=NULL
for (cell in cell_types) {
  for (n in which(cell_types==cell):length(cell_types)) {
    from <- c(from, cell_types[n])
    to<-c(to, cell)
  }
}
perm=data.frame(from=from, to=to)

res <- 
  perm %>%
  purrr::pmap(
    function(from, to){
      iTALK::FindLR(
        cell_res[[from]],
        cell_res[[to]],
        datatype = "DEG",
        comm_type = "intercell",
        database = db_italk)
    }
  ) %>% bind_rows() 

res$Pair.Name=paste(res$ligand,res$receptor, sep="_")

# Interaction funtional annotations ####
res= left_join(res, db_italk %>% dplyr::select(c(1:6,12:14,21,22,48))) %>% filter(curation_effort > 0)

# Intercell annotation https://r.omnipathdb.org/articles/omnipath_intro.html#intercell
intercell <- import_omnipath_intercell(scope = "generic",
                                       proteins = unique(c(res$ligand, res$receptor)))

intercell_ligand= intercell[match(res$ligand,intercell$genesymbol),] %>% 
  dplyr::select(1,2) %>% dplyr::rename_with(~ paste("ligand", .x, sep = "_"))

intercell_target= intercell[match(res$receptor,intercell$genesymbol),] %>% 
  dplyr::select(1,2) %>% dplyr::rename_with(~ paste("target", .x, sep = "_"))

res=bind_cols(res, intercell_ligand,intercell_target) %>% 
  relocate(category_intercell_source,category_intercell_target,ligand_category, target_category, .after = Pair.Name)

path_annot=c("SignaLink_pathway","KEGG-PC","CellChatDB")

for (annotation in path_annot){
  res=anot_path(res, annotation)
}

# ####

res <- res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]
# circos=rbind(circos,res)
if(is.null(res)){
  print('No significant pairs found')
}else{
  circos=res
  circos=rbind(circos,res)
  if(nrow(res)>=20){circos=circos[1:25,]}
  gene_col<-structure(viridis::plasma(length(unique(c(circos$ligand_category,circos$target_category))), direction = 1),names=unique(c(circos$ligand_category,circos$target_category)))
  gene_col <- structure(c(gene_col[circos$ligand_category], gene_col[circos$target_category]), 
                        names = c(circos$ligand, circos$receptor))
  
  jpeg(paste0(tissue, "/Comunicacion_celular/0logfc_SDID_CCC.jpeg"),width=7, height=6, units="in",res=250)
  par(mfrow=c(1,2))
  LRPlot(circos,datatype='DEG',gene_col=gene_col,cell_col=cell_col,
         link.arr.lwd=circos$cell_from_logFC,link.arr.width=circos$cell_to_logFC)
  NetView(res,col=cell_col[unique(c(arrange(res, cell_from)$cell_from,arrange(res, cell_to)$cell_to))],
          vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  title(paste("SDID differential interactions in",ts), cex.main= .8)
  dev.off()
}

write.table(circos,paste0(tissue, "/Comunicacion_celular/0logfc_SDID_CCC_intercell.tsv"), quote = F, col.names = T, row.names = F, sep="\t")
