# Single cell RNA-seq data analysis
# Plot marker expression
# Autor: Fernando Gordillo Gonzalez
# Date: 31/07/2023
# Args: nombre_estudio ncols_plot marker_genes

library(SingleCellExperiment)
library(celda)
library(scater)

#---------------------- Make directories & Load data --------------------------####

setwd("/home/fgordillo/sc-RNAseq")
# Set wd
args=commandArgs(trailingOnly = TRUE)
estudio = as.character(args[1])
ncol = as.numeric(args[2])

# Create directories
setwd(paste0(estudio,"/figures/5.Cell_annotation/"))
getwd()

# Load data
sce = readRDS(paste0("../../results/5.Cell_annotation/sce_",tolower(estudio),"_annot_filt.rds"))

# plotExpression.local function  ####
plotExpression.local <- function (object, features, x = NULL, exprs_values = "logcounts", 
                                  log2_values = FALSE, colour_by = NULL, shape_by = NULL, 
                                  size_by = NULL, by_exprs_values = exprs_values, xlab = NULL, 
                                  feature_colours = TRUE, one_facet = TRUE, ncol = 2, scales = "fixed", 
                                  other_fields = list(), swap_rownames = NULL, ...) 
{
  if (!is(object, "SingleCellExperiment")) {
    stop("object must be an SingleCellExperiment object.")
  }
  if (exprs_values == "exprs" && !(exprs_values %in% assayNames(object))) {
    exprs_values <- "logcounts"
  }
  exprs_vals <- vector("list", length(features))
  for (i in seq_along(features)) {
    current <- retrieveCellInfo(object, features[i], search = c("assays", 
                                                                "altExps"), exprs_values = exprs_values, swap_rownames = swap_rownames)
    features[[i]] <- current$name
    if (is.null(current$value)) {
      stop("cannot find '%s' in 'object'", features[i])
    }
    exprs_vals[[i]] <- unname(current$value)
  }
  nfeatures <- length(features)
  if (log2_values) {
    exprs_val <- lapply(exprs_vals, function(x) log2(x + 
                                                       1))
    ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
  }
  else {
    ylab <- paste0("Expression (", exprs_values, ")")
  }
  evals_long <- data.frame(Feature = rep(factor(features, 
                                                features), lengths(exprs_vals)), Y = unlist(exprs_vals))
  x_by_out <- retrieveCellInfo(object, x, exprs_values = exprs_values)
  xcoord <- x_by_out$val
  if (is.null(xlab)) {
    xlab <- x_by_out$name
  }
  evals_long$X <- rep(xcoord, nfeatures)
  vis_out <- scater:::.incorporate_common_vis_col(evals_long, se = object, 
                                                  colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                                  by_exprs_values = by_exprs_values, other_fields = other_fields, 
                                                  multiplier = rep(seq_len(ncol(object)), nfeatures), 
                                                  swap_rownames = swap_rownames)
  evals_long <- vis_out$df
  colour_by <- vis_out$colour_by
  shape_by <- vis_out$shape_by
  size_by <- vis_out$size_by
  if (is.null(evals_long$X)) {
    evals_long$X <- evals_long$Feature
  }
  else {
    one_facet <- FALSE
  }
  feature_colours <- (feature_colours && one_facet && is.null(colour_by))
  if (feature_colours) {
    evals_long$fill_by <- evals_long$colour_by <- evals_long$Feature
    fill_by <- colour_by <- "Feature"
  }
  else {
    fill_by <- NULL
  }
  plot_out <- scater:::.central_plotter(evals_long, xlab = xlab, ylab = ylab, 
                                        shape_by = shape_by, colour_by = colour_by, size_by = size_by, 
                                        fill_by = fill_by, ..., point_FUN = NULL)
  if (!one_facet) {
    plot_out <- plot_out + facet_wrap(~Feature, ncol = ncol, 
                                      scales = scales, strip.position="right")
  }
  if (is.null(x)) {
    plot_out <- plot_out + theme(axis.text.x = element_text(angle = 60, 
                                                            vjust = 1, hjust = 1), axis.ticks.x = element_blank(), 
                                 plot.margin = unit(c(0.03, 0.02, 0.05, 0.02), "npc"))
  }
  if (feature_colours) {
    plot_out <- plot_out + guides(fill = "none", colour = "none")
  }
  plot_out
}

#---------------------------------- Plot -----------------------------------####

# Paper markers
ref_markers=strsplit(args[3], ",")[[1]]
ref_markers=ref_markers[which(ref_markers %in% rowData(sce)$SYMBOL)]
rownames(sce)=rowData(sce)$SYMBOL

jpeg(file="after_filt/Markers_violin_plots.jpg", width=6*ncol, height=15, units="in", res=300)
plotExpression(sce, features = ref_markers, x = "cell_type", colour_by = "cell_type", ncol = ncol) + xlab("Cell types") + ylab("Normalised gene expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top",
        strip.text.x = element_text(size = 13)) +
  scale_color_viridis_d(option = "H")
dev.off()


jpeg(file="after_filt/Markers_violin_plots_clusters.jpg", width=8, height=2*length(ref_markers), units="in", res=300)
plotExpression.local(sce, features = ref_markers, x = "cluster", colour_by = "cell_type", ncol = 1) + xlab("Cluster") + ylab("Normalised gene expression") +
  theme(legend.position = "top",
        strip.text.y = element_text(size = 13),
        strip.background = element_rect(color="black", fill="#e5e5e5", size=0.5)) +
  scale_color_viridis_d("Cell type",option = "H")
dev.off()


# plotGroupedHeatmap(sce, features=ref_markers, group="cluster", 
#                    center=TRUE, zlim=c(-3, 3),swap_rownames="SYMBOL",cluster_rows=F)

a=ceiling(sqrt(length(ref_markers)))
jpeg(file="after_filt/Markers_UMAP_plots.jpg", width=5*a, height=5*a, units="in", res=250)
umap <- reducedDim(sce, "UMAP_300")
plotDimReduceFeature(logcounts(sce),
                     dim1 = umap[, 1],
                     dim2 = umap[, 2],
                     features = c(ref_markers),
                     ncol = a)
dev.off()

