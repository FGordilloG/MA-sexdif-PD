# MA-sexdif-PD

Meta-analysis of sex-specific differential expression in single-nucleus RNA-seq datasets of post-mortem substantia nigra pars compacta (SN) and cortical (CT) samples from Parkinson's disease (PD) patients and controls.

This repository contains the full analysis code used in the paper, in the order in which it was executed. Results can be explored interactively in the **PATOSS** web application ([link]).

- **Paper:** [citation / DOI]
- **Archived release:** [Zenodo DOI]
- **Data:** all datasets are public. Accession numbers, sample sizes and donor metadata are listed in Table 1 of the paper.

---

## How the analysis was run

Scripts are numbered in execution order. Each script is run from the command line and takes its arguments in the order given below, for example:

```bash
Rscript 1.QC.R GSE157783 <ncells> <UMIs> <ngenes> <%mito>
```

The parameter values used for every dataset (quality control thresholds, clustering parameters, number of principal components and model covariates) are reported in Supplementary Table [SX] of the paper.

### Stage 1 — Individual analysis (run once per dataset)

| # | Script | What it does | Arguments |
|---|---|---|---|
| 1 | `1.QC.R` | Quality control: removal of empty droplets, low-quality nuclei, doublets and lowly detected genes | `study_ID ncells UMIs ngenes %mito` |
| 2 | `2.Normalization.R` | Deconvolution-based size factor normalisation and log-transformation | `study_ID Ensembl_ID_plot` |
| 3 | `3.Feature_selection.R` | Selection of the most variable genes | `study_ID` |
| 4 | `4.Dim_reduct.R` | PCA, selection of the number of components and UMAP. Optionally applies batch correction (see below) | `study_ID batch_correction covariates` |
| 5 | `5.Clustering.R` | Cell clustering combining k-means and shared nearest neighbours | `study_ID centers k` |
| 6 | `6.Cell_type_annot.R` | Cell type annotation by consensus of four methods (SingleR, CHETAH, SCINA, BRETIGEA) | `study_ID sce_ref_path genelist_path recap` |
| 7 | `7.Plot_marker_expression.R` | Verification of the annotation through the expression of canonical marker genes | `study_ID ncols_plot marker_genes` |

### Stage 2 — Differential analysis (run per dataset and cell type)

Two branches run in parallel on the annotated objects produced in stage 1.

| # | Script | What it does | Arguments |
|---|---|---|---|
| 8 | `8.DE_MAST.R` | Gene-level differential expression with MAST, using the condition-by-sex interaction (SDID contrast) | `study_ID cell_type gene_lengths_path` |
| 13 | `13.hipathia.R` | Signalling pathway analysis: computation of sub-pathway activation values per nucleus | `study_ID` |
| 14 | `14.hipathia_DEA.R` | Differential activation of sub-pathways with MAST, using the same contrast | `study_ID cell_type` |

No inference is drawn from the individual datasets at this stage. Each analysis provides an effect size (logFC) and its standard error, which are the input of the meta-analysis.

### Stage 3 — Integration

| # | Script | What it does | Arguments |
|---|---|---|---|
| 9 | `9.Meta-analysis.R` | Random-effects meta-analysis (DerSimonian–Laird) of the gene-level and pathway-level results, per cell type and region, with Benjamini–Hochberg correction and diagnostic tests | see header of the script |

### Stage 4 — Downstream analyses (run on the meta-analysis results)

| # | Script | What it does | Arguments |
|---|---|---|---|
| 10 | `10.Functional_enrichment.R` | GO Biological Process enrichment with the weight01 algorithm (topGO) | `cell_type tissue` |
| 11 | `11.Semantic_GO_clustering.R` | Semantic clustering of the enriched GO terms and summary heatmaps | `cell_type tissue` |
| 12 | `12.CCC_italk.R` | Cell–cell communication analysis (iTALK) using ligand–receptor pairs from OmniPath | `cell_type contrast` |

---

## Batch correction

**Batch correction was not applied in the analyses reported in the paper.**

Each dataset was processed and analysed separately, and the variability between studies is accounted for by the random-effects meta-analysis, so no correction across datasets is required by the design.

Within datasets, batch correction was nevertheless tested during the development of the pipeline, using the Harmony implementation provided in `Batch_correction.R`. It did not produce relevant changes in the resulting embeddings, clustering or cell type annotation, and it was therefore left out of the final pipeline.

The function is kept in the repository, and is still reachable through the `batch_correction` argument of `4.Dim_reduct.R`, for two reasons: so that the code can be reused in settings where correction is needed, and so that this decision can be re-examined by anyone wishing to reproduce or extend the analysis.

---

## Requirements

Analyses were run in R [version]. The version of every package used is listed in Supplementary Table [SY] of the paper and in `sessionInfo.txt`.

Main packages: SingleCellExperiment, DropletUtils, scater, scran, bluster, scDblFinder, uwot, PCAtools, SingleR, CHETAH, SCINA, BRETIGEA, MAST, metafor, hipathia, topGO, GOSemSim, simplifyEnrichment, STRINGdb, iTALK, OmnipathR, harmony, shiny.

## Contact

Fernando Gordillo González — [email / institution]
