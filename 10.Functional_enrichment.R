# sn-RNA-seq sex-difs atlas in PD
# Functional profiling with weight01
# Autor: Fernando Gordillo Gonzalez
# Args: cell_type tissue

#----------------------------- Load packages ----------------------------------####

library(GOSemSim)
library(genefilter)
library(GEOquery)
library(hgu133plus2.db)
library(topGO)
library(Hmisc)
library(limma)
library(dplyr)
library(GO.db)
library(simplifyEnrichment)
library(rrvgo)
library(STRINGdb)

#---------------------- Make directories & Load data --------------------------####

setwd("~/PDsex_atlas/Integracion/MA")
getwd()

# Arguments
args=commandArgs(trailingOnly = TRUE)
cell=as.character(args[1])
tissue=as.character(args[2])

# Load data
dataset = read.delim(paste("../",tissue,cell,"FM/all.genes.tsv",sep="/"), row.names = 1)

## Global variables
# Generic words to avoid in the word-clouds
no_words=c("regulation", "process","cell", "signaling", "pathway","particle","cellular", 
           "nuclear","negative", "positive","tissue", "activity","cardiac","muscle")
medida = "Rel"
organismo = "org.Hs.eg.db"
algoritmo = "personalizado"
pvalor=0.05

# Init db for string PPI generations
string_db <- STRINGdb$new( version="11.5", species=9606,
                           score_threshold=400) 

# Global GO variables for enrichment and redundance prunning
GO_anc = as.list(GOBPANCESTOR)
GO_def = suppressMessages(AnnotationDbi::select(GO.db,
                                                keys=(keys(GO.db)),
                                                columns=c("GOID","TERM","DEFINITION"),
                                                keytype="GOID"))
rownames(GO_def) = GO_def$GOID
GO_semsim = godata(org.Hs.eg.db, ont="BP")


#---------------------- Functions --------------------------####


filtrado_genes <- function(dataset, pvalor =  0.05, fc = 0.5) {
  
  # Filtrar genes
  sig = dataset[dataset$p.adjust.fdr < pvalor,]
  down = rownames(sig[(sig$logFC < (-fc)),])
  up = rownames(sig[(sig$logFC > (fc)),])
  geneList = as.vector(dataset$p.adjust.fdr)
  names(geneList) = row.names(dataset)
  
  # Escogemos los genes up y hacemos análisis
  upList <- factor(as.integer(names(geneList) %in% up))
  names(upList) <- names(geneList)
  
  # Escogemos los genes down y hacemos analisis
  downList <- factor(as.integer(names(geneList) %in% down))
  names(downList) <- names(geneList)
  
  print(table(upList))
  print(table(downList))
  
  return(list(upList, downList))
}

weight_enrich <- function(genes_sig,dataset, pvalor =  0.05, organismo = "org.Hs.eg.db", ontology = "BP", notacion = "ENSEMBL"){
  pre_enriquecimiento <- new("topGOdata", 
                                ontology = ontology, 
                                allGenes = genes_sig, 
                                geneSel = function(x)(x == 1), 
                                nodeSize = 15,
                                annot = annFUN.org,
                                mapping = organismo,
                                ID = notacion)
  
  ## 3.2.2 Perform the analysis
  weight01_up <- runTest(pre_enriquecimiento, algorithm = "weight01", statistic = "fisher", cutOff = pvalor)
  
  weight01_t_up = GenTable(pre_enriquecimiento, weight01 = weight01_up,  topNodes = length(score(weight01_up)))
  weight01_t_up$TERM = GO_def[weight01_t_up$GO.ID, "TERM"]
  weight01_t_up$p.value = sort(weight01_up@score) # Esto lo hacemos para incorporar los números más decimales (weight01 corresponde con los p.valores)
  weight01_t_up$weight01=NULL
  
  ## 3.2.3. Add gene info to the data.frame
  sigGenes <- sigGenes(pre_enriquecimiento)
  AnnoList <- lapply(weight01_t_up$"GO.ID", function(x) as.character(unlist(genesInTerm(object = pre_enriquecimiento, whichGO = x))))
  SigList <- lapply(AnnoList, function(x) intersect(x, sigGenes))
  SigList <- lapply(SigList, function(x) dataset[rownames(dataset) %in% x,"SYMBOL"])
  weight01_t_up$"SYMBOLS" <- sapply(SigList, paste, collapse = ",")
  
  # To check propagation
  # print(showSigOfNodes(pre_enriquecimiento, weight01_up@score, firstSigNodes = 20, useInfo = "all"))
  
  return(weight01_t_up)
}

register_clustering_methods(
  personalizado = function(mat) binary_cut(mat, try_all_partition_fun = TRUE, cutoff = 0.8)
)

postprocesado <- function(a,tissue,cell,label,  medida = "Rel", algoritmo = "personalizado", umbral = 0.7, organismo = 'org.Hs.eg.db') {
  go_id = names(a)
  
  # Calculamos reducción con simplifyEnrichment
  mat = GO_similarity(names(a), ont = "BP", db = organismo , measure = medida)
  df_fisher = simplifyGO(mat, method = algoritmo, plot = FALSE)
  jpeg(paste0("Enrichment_SE_clust_",cell,label,".jpg"), width = 9, height = 6, units = "in", res=250)
    ht_clusters(mat, df_fisher$cluster, draw_word_cloud = T, order_by_size = T, min_term = 3, exclude_words=no_words)
  dev.off()
  
  
  # Calculamos reducción con Revigo
  simMatrix = calculateSimMatrix(names(a), orgdb = organismo, ont = "BP", method = medida)
  scores <- setNames(-log10(a), names(a)) #REDUCIR A -log10 PORQUE CUANTO MAS GRANDE MEJOR
  reducedTerms = reduceSimMatrix(simMatrix, scores, threshold = umbral, orgdb = organismo)
  jpeg(paste0("Enrichment_RVGO_tree_",cell,label,".jpg",sep=""), width = 12, height = 10, units = "in", res=250)
    treemapPlot(reducedTerms)
  dev.off()
  jpeg(paste0("Enrichment_RVGO_clust_",cell,label,".jpg"), width = 9, height = 6, units = "in", res=250)
    ht_clusters(simMatrix, reducedTerms$cluster, draw_word_cloud = T, order_by_size = T, min_term = 3, column_title = paste("Enrichment RVGO",tissue,cell, label), exclude_words=no_words)
  dev.off()
  
  output_revigo = reducedTerms[,c(1,7,2)]
  colnames(output_revigo) = colnames(df_fisher)
  # Creamos un dataframe donde asignamos a cada término nuestro todos sus ancestros
  x<-stack(GO_anc[go_id])
  colnames(x) = c("ancestral", "go_id")
  
  # Combinamos nuestro resultado de simplifyEnrichment con nuestro dataframe de ancestros para tener la columna "clusters"
  x2<-left_join( x, df_fisher[,c(3,1)], by = setNames("id", "go_id"))
  x3<-left_join( x, output_revigo[,c(3,1)], by = setNames("id", "go_id"))
  x3 = x3[complete.cases(x3),]
  
  # Agrupamos creando una dupla entre los ancestros y el cluster en el que aparece de forma que (sigue abajo)
  x2<-group_by(x2,cluster, ancestral)
  x3<-group_by(x3,cluster, ancestral)
  cluster_size = (as.numeric(as.matrix(table(add_count(group_by(df_fisher,cluster))$cluster))))
  cluster_size_revigo = (as.numeric(as.matrix(table(add_count(group_by(output_revigo,cluster))$cluster))))
  
  # Contamos las veces que un ancestro aparece en un cluster concreto, es como responder a la pregunta:
  # ¿Cuántas veces aparece el ancestro X en el cluster Y?, la respuesta nos la da la nueva columna "n"
  x2<-add_count(x2) 
  x3<-add_count(x3) 
  x2<-ungroup(x2)
  x3<-ungroup(x3)
  
  # Creamos nuevas variables para nuestro dataframe
  # Usamos la variable GO_semsim que contiene todos los datos de IC del árbol GO par BP
  x2$IC = GO_semsim@IC[x2$ancestral]
  x3$IC = GO_semsim@IC[x3$ancestral]
  
  # Creamos de nuevo la columna n que nos dice cuantas veces sale cada término padre en un cluster
  x2<-group_by(x2,cluster, ancestral)
  x3<-group_by(x3,cluster, ancestral)
  
  # Calculamos la importancia de cada término padre en función de las veces que sale y su IC
  x2$importance = ((x2$n/cluster_size[x2$cluster])*x2$IC^2)
  x3$importance = ((x3$n/cluster_size_revigo[x3$cluster])*x3$IC^2)
  
  # Seleccionamos para cada cluster el mejor candidato
  x2 = x2 %>% group_by(cluster) %>% dplyr::slice(which.max(importance))
  x3 = x3 %>% group_by(cluster) %>% dplyr::slice(which.max(importance))
  
  # Añadimos informacion de los padres seleccionados con el objeto "GO_def" cargado al principio
  x2 <- left_join(x2, GO_def[,c(1,2)], by = setNames("GOID", "ancestral"))
  x3 <- left_join(x3, GO_def[,c(1,2)], by = setNames("GOID", "ancestral"))
  x2 <- x2[,c(1,2,3,4,5,6,7)]
  x3 <- x3[,c(1,2,3,4,5,6,7)]    
  x2$cluster_size = cluster_size
  x3$cluster_size = cluster_size_revigo
  
  for(i in 1:nrow(x2)) { 
    s <- x2[i,8]
    if(s$cluster_size==1) {
      patata = left_join(x2[i,], GO_def[,c(1,2)], by = setNames("GOID", "go_id"))
      patata[,7]=patata[,9]
      x2[i,] = patata[,-9]
      #x2[i,5] = (length(which(x2$go_id==x2[i,2]$go_id)))*(GO_semsim@IC[x2[i,2]$go_id]^2)
      x2[i,1]=x2[i,2]
      #x2[i,4]=as.numeric(GO_semsim@IC[as.character(x2[i,2])])
    }
  }
  
  for(i in 1:nrow(x3)) { 
    s <- x3[i,8]
    if(s$cluster_size==1) {
      patata = left_join(x3[i,], GO_def[,c(1,2)], by = setNames("GOID", "go_id"))
      patata[,7]=patata[,9]
      x3[i,] = patata[,-9]
      #x3[i,5] = (length(which(x3$go_id==x3[i,2]$go_id)))*(GO_semsim@IC[x3[i,2]$go_id]^2)
      x3[i,1]=x3[i,2]
      #x3[i,4]=as.numeric(GO_semsim@IC[as.character(x3[i,2])])
    }
  }
  
  x2$n =   (x2$n/x2$cluster_size)*100
  x3$n =   (x3$n/x3$cluster_size)*100
  x2 = x2[,-4]
  x3 = x3[,-4]
  
  # Creamos la combinacion de ambos
  combinacion<-merge(x2, x3, by = intersect("TERM", "TERM"))
  combinacion<-unique(combinacion[,-c(3,4,5,8,9,10,11)])
  colnames(combinacion)=c("TERM","Ancestro","Importanica-SE","cluster_size_SE","Importancia-R","cluster_size_revigo")
  
  return(list(simplifyenrichment = df_fisher, revigo = reducedTerms, clusterizacion_se = x2, clusterizacion_revigo = x3, combinacion = combinacion))
  
}


SE <- function(weight01_t_up, cell, tissue, label, pvalor =  0.05, organismo = "org.Hs.eg.db", medida="Rel", no_words=character(0)){
  mat = GO_similarity(weight01_t_up$GO.ID, ont = "BP", db = organismo , measure = medida)
  df_fisher = simplifyGO(mat, method = "personalizado", plot = F)
  jpeg(paste0("Enrichment_SE_clust_",cell,label,".jpg"), width = 9, height = 6, units = "in", res=250)
    ht_clusters(mat, df_fisher$cluster, draw_word_cloud = T, order_by_size = T, min_term = 3, column_title = paste("Enrichment",tissue,cell, label), exclude_words=no_words)
  dev.off()
}

revigo <- function(dt,label, pvalor =  0.05, organismo = "org.Hs.eg.db", medida="Rel"){
  simMatrix = calculateSimMatrix(dt$GO.ID, orgdb = organismo, ont = "BP", method = medida)
  scores <- setNames(-log10(dt$p.value), dt$GO.ID)
  reducedTerms = reduceSimMatrix(simMatrix, scores, orgdb = organismo)
  
  jpeg(paste0("Enrichment_RVGO_tree_",cell,label,".jpg",sep=""), width = 12, height = 10, units = "in", res=250)
  treemapPlot(reducedTerms)
  dev.off()
}

my_string<-function(res, label,cutoff=0.05, notation="ENSEMBL"){
  res_string <- string_db$map( res, notation, removeUnmappedRows = TRUE )
  # plot_network
  string_db$get_png( res_string$STRING_id, file=paste0("STRING_",cell,label,".png"))
  
  
  # enrichment
  enrichment <- string_db$get_enrichment( res_string$STRING_id )
  
  if(nrow(enrichment)>0){
    write.table(enrichment[,c(9,10,1:4,8)], file=paste0("STRING_enrich_",cell,label,".tsv"), sep="\t", quote = TRUE, row.names = F)
  }

  # get url
  url <- string_db$get_link(res_string$STRING_id)
  cat(paste(paste0(cell,label),url, sep="\t"),"\n" , file = paste0("../STRING_urls.tsv") ,append = T)
}

write_res <- function(weight01_t_up, res, cell, label){
  weight01_t_up$SE_cluster <- res[["simplifyenrichment"]][match(weight01_t_up$GO.ID, res[["simplifyenrichment"]]$id),"cluster"]
  weight01_t_up$RVGO_cluster <- res[["revigo"]][match(weight01_t_up$GO.ID, res[["revigo"]]$go),"cluster"]
  write.table(weight01_t_up, file=paste0("Enrichment_weight01_",cell,label,".tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  
  cat("simplifyEnrichment\n", file = paste0("Enrichment_abstract_",cell,label,".tsv"))
  write.table(res[["clusterizacion_se"]] %>% arrange(desc(cluster_size)), file=paste0("Enrichment_abstract_",cell,label,".tsv"), append = TRUE, col.names=TRUE, row.names=FALSE, sep="\t")
  cat("\n\nREVIGO\n", file = paste0("Enrichment_abstract_",cell,label,".tsv"), append = TRUE)
  write.table(res[["clusterizacion_revigo"]] %>% arrange(desc(cluster_size)), file=paste0("Enrichment_abstract_",cell,label,".tsv"), append = TRUE, col.names=TRUE, row.names=FALSE, sep="\t")
  cat("\n\nCombinacion\n", file = paste0("Enrichment_abstract_",cell,label,".tsv"), append = TRUE)
  write.table(res[["combinacion"]], file=paste0("Enrichment_abstract_",cell,label,".tsv"), append = TRUE, col.names=TRUE, row.names=FALSE, sep="\t")
}


# ------------------------ Enrichment ----------------------------####

numbers=matrix(NA, nrow = length(cell_types), ncol = length(tissue)*2)
rownames(numbers)=cell_types
colnames(numbers) <- as.character(rbind(paste0(tissue,"_up"),paste0(tissue,"_down")))


# Select significant genes: p.adj < 0.05 and |logFC| > 0.5
genes_sig_up <- filtrado_genes(dataset = dataset)[[1]] # Genes up
genes_sig_down <- filtrado_genes(dataset = dataset)[[2]]  #Genes down

dir.create(paste(tissue,cell,sep="/"))
setwd(paste(tissue,cell,sep="/"))

if (sum(genes_sig_up == 1) > 0) {
cat(names(genes_sig_up[genes_sig_up == 1]),sep="\n", file = paste0("Genes_",cell,"_up.txt"))
  
  # weight01 enrich
  weight01_t_up=weight_enrich(genes_sig_up,dataset)

  # weight01_t_up=read.delim(paste0("Enrichment_weight01_",cell,"_up.tsv"))
  
  # Clustering GOs
  if (sum(weight01_t_up$p.value<pvalor) > 5) {
    weight01 = weight01_t_up[weight01_t_up$p.value<pvalor,]
    weight01 = structure(weight01$p.value, names=weight01$GO.ID)
    res=postprocesado(weight01, tissue, cell, "_up")
    write_res(weight01_t_up, res, cell, "_up")
  } else {
    write.table(weight01_t_up, file=paste0("Enrichment_weight01_",cell,"_up.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  }
  numbers[cell, paste0(tissue,"_up")] = sum(weight01_t_up$p.value<pvalor)
  # STRING
  my_string(data.frame(ENSEMBL=names(genes_sig_up[genes_sig_up == 1])), "_up")
}

if (sum(genes_sig_down == 1) > 0) {
cat(names(genes_sig_down[genes_sig_down == 1]),sep="\n", file = paste0("Genes_",cell,"_down.txt"))
  
  # weight01 enrich
  weight01_t_down=weight_enrich(genes_sig_down,dataset,"_down")
  # weight01_t_down=read.delim(paste0("Enrichment_weight01_",cell,"_down.tsv"))
  
  # Clustering GOs
  if (sum(weight01_t_down$p.value<pvalor) > 5) {
    weight01 = weight01_t_down[weight01_t_down$p.value<pvalor,]
    weight01 = structure(weight01$p.value, names=weight01$GO.ID)
    res=postprocesado(weight01, tissue, cell, "_down")
    write_res(weight01_t_down, res, cell, "_down")
  } else {
    write.table(weight01_t_down, file=paste0("Enrichment_weight01_",cell,"_down.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  }
  numbers[cell, paste0(tissue,"_down")] = sum(weight01_t_down$p.value<pvalor)
  # STRING
  my_string(data.frame(ENSEMBL=names(genes_sig_down[genes_sig_down == 1])), "_down")
}
