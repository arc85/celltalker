## Cluster and filter HCA BM dataset from Seurat Data
## Mar 28 2022

## Load packages
library(celltalker)
library(Seurat)
library(SeuratData)
library(tidyverse)

## Load CMBC data
data(hcabm40k)

## Data processing
cbmc <- hcabm40k %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(cbmc)

cbmc <- RunUMAP(cbmc,dims=1:8) %>%
  FindNeighbors(.,dims=1:8) %>%
  FindClusters(.,res=0.7)

## ID cell types
# Check out clusters
DimPlot(cbmc,group.by="RNA_snn_res.0.7",label=T)

FeaturePlot(cbmc,c("CD34","CD14","FCGR3A","MS4A1","CD1C","IL3RA"))
FeaturePlot(cbmc,c("CD3D","CD4","CD8A","FCGR3A","HBB"))

# Assign cell types
cell_types <- vector("logical",length=ncol(cbmc))
names(cell_types) <- colnames(cbmc)

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="0" |
  cbmc@meta.data$RNA_snn_res.0.7=="2" |
  cbmc@meta.data$RNA_snn_res.0.7=="6" ] <- "CD4 T cells"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="3"] <- "CD8 T cells"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="5"] <- "NK cells"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="8" |
  cbmc@meta.data$RNA_snn_res.0.7=="11"] <- "RBCs"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="1"] <- "CD14 monocytes"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="14"] <- "CD16 monocytes"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="15"] <- "CD1C DC"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="10"] <- "CD34 HPSC"

cell_types[cbmc@meta.data$RNA_snn_res.0.7=="4"] <- "B cells"

cbmc[["cell_types"]] <- cell_types

# Filter cells - remove transitional cell populations and RBCs
cbmc_filtered <- cbmc[,!cbmc[["cell_types"]]=="FALSE"]
cbmc_filtered <- cbmc_filtered[,!cbmc_filtered[["cell_types"]]=="RBCs"]

# Check out final object
DimPlot(cbmc_filtered,group.by="cell_types")

# Save UMAP embeddings and cell types
hca_bm_umap_cell_types <- data.frame(Embeddings(cbmc_filtered,reduction="umap"),
  cell_types=cbmc_filtered[["cell_types"]])

## Save data
usethis::use_data(hca_bm_umap_cell_types)
