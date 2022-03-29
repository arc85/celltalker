
<!-- README.md is generated from README.Rmd. Please edit that file -->

# celltalker

The goal of celltalker is to infer putative ligand and receptor
interactions from single-cell RNAseq data. This is accomplished by
evaluating known cognate ligand/receptor interactions across groups of
cells. Interactions are scored by jointly weighting the expression
levels of ligands and receptors, and significance is evaluated by
comparing to a background distribution of scrambled ligands and
receptors.

A recent refactoring of this package has made the interface much simpler
to use. We also provide a wrapper to the functionality of the [circlize
package](https://jokergoo.github.io/circlize_book/book/) for creating
circos plots of ligands and receptors.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("arc85/celltalker")
```

## TL;DR example use

Here’s a basic use case for celltalker based on a subset of the 40,000
bone marrow single-cell dataset available from the Human Cell Atlas.
Check out SeuratData [here](https://github.com/satijalab/seurat-data) to
download the example dataset.

Note that we have identified canonical immune cell types and filtered
the dataset for ease of analysis. A vignette documenting this processing
is coming soon, but please utilize the **hca_bm_umap_cell_types**
data.frame for the vignette below.

``` r
suppressMessages({
library(celltalker)
library(Seurat)
library(SeuratData)
library(tidyverse)
})
#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

#> Warning in if (is.na(desc)) {: the condition has length > 1 and only the first
#> element will be used

data(hcabm40k)

# Filter cell and assign cell types to the dataset
hca_bm <- hcabm40k[,rownames(hca_bm_umap_cell_types)]
hca_bm[["cell_types"]] <- hca_bm_umap_cell_types$cell_types

# Process data
hca_bm <- NormalizeData(hca_bm)

# Add UMAP coordinates
hca_bm[["umap"]] <- CreateDimReducObject(embeddings=as.matrix(hca_bm_umap_cell_types[,1:2]),key="UMAP_")
#> Warning: No assay specified, setting assay as RNA by default.

# View cell types
DimPlot(hca_bm,group.by="cell_types")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
## Run celltalker
hca_bm_interactions <- celltalk(input_object=hca_bm,
  metadata_grouping="cell_types",
  ligand_receptor_pairs=ramilowski_pairs,
  number_cells_required=100,
  min_expression=1000,
  max_expression=20000,
  scramble_times=10)

## Identify top statistically significant interactions
top_stats <- hca_bm_interactions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05) %>%
  group_by(cell_type1) %>%
  top_n(3,interact_ratio) %>%
  ungroup()

## Generate a circos plot
colours_use <- RColorBrewer::brewer.pal(n=length(unique(hca_bm$cell_types)),"Set2")

suppressMessages(
circos_plot(ligand.receptor.frame=top_stats,colors=colours_use,lig.col="blue",rec.col="red")
)
```

<img src="man/figures/README-example-2.png" width="100%" />
