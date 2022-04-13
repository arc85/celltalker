## Preprocessing of samples for celltalker analysis
## Dec 22 2021

## Load Seurat
library(Seurat)

## Read in samples and create metdata
files.to.read <- grep("\\.R",list.files("data-raw"),invert=TRUE,value=TRUE)
ser.list <- vector("list",length=length(files.to.read))

for (i in 1:length(files.to.read)) {

  tmp.sample <- Read10X(paste("data-raw",files.to.read[i],sep="/"))
  tmp.meta <- data.frame(sample_id=rep(NA,ncol(tmp.sample)),sample_type=rep(NA,ncol(tmp.sample)))
  tmp.meta$sample_id <- rep(files.to.read[i],length(tmp.meta$sample_id))
  tmp.meta$sample_type <- ifelse(grepl("PBMC",tmp.meta$sample_id),"PBMC","Tonsil")
  rownames(tmp.meta) <- colnames(tmp.sample)

  ser.list[[i]] <- CreateSeuratObject(tmp.sample,meta.data=tmp.meta)

}

## Merge into one Seurat Object
ser <- merge(ser.list[[1]],ser.list[2:length(ser.list)])

## Seruat clustering and visualization workflow
ser <- NormalizeData(ser)
ser <- FindVariableFeatures(ser)
ser <- ScaleData(ser)
ser <- RunPCA(ser)
ElbowPlot(ser)

ser <- RunUMAP(ser,dims=1:10)
ser <- FindNeighbors(ser,dims=1:10)
ser <- FindClusters(ser,res=0.3)

## Identify cell types
DimPlot(ser,group.by="sample_type")
DimPlot(ser,label=T)
FeaturePlot(ser,c("CD3D","CD8A","CD4",
                  "CD14","FCGR3A",
                  "MS4A1","MZB1",
                  "IL3RA","CLEC10A",
                  "HBB","PPBP"))

## Add cell type metadata
data.to.add <- vector("logical",length=ncol(ser))

data.to.add[ser@meta.data$RNA_snn_res.0.3=="0" |
              ser@meta.data$RNA_snn_res.0.3=="6" |
              ser@meta.data$RNA_snn_res.0.3=="7"] <- "B cells"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="10"] <- "Plasmablasts"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="1" |
              ser@meta.data$RNA_snn_res.0.3=="2" |
              ser@meta.data$RNA_snn_res.0.3=="3" |
              ser@meta.data$RNA_snn_res.0.3=="5" |
              ser@meta.data$RNA_snn_res.0.3=="8" ] <- "NK and T cells"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="4"] <- "CD14+CD16- monocytes"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="9"] <- "CD14-CD16+ monocytes"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="11"] <- "CD1C+ DCs"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="12"] <- "pDCs"

data.to.add[ser@meta.data$RNA_snn_res.0.3=="13"] <- "RBCs"

ser[["cell_types"]] <- data.to.add

## Visualize cell types
DimPlot(ser,group.by="cell_types",split.by="sample_type")

## Subcluster to identify NK and T cells
ser.t <- ser[,ser@meta.data$cell_types=="NK and T cells"]
ser.t <- FindVariableFeatures(ser.t)
ser.t <- ScaleData(ser.t)
ser.t <- RunPCA(ser.t)
ElbowPlot(ser.t)

ser.t <- RunUMAP(ser.t,dims=1:10)
ser.t <- FindNeighbors(ser.t,dims=1:10)
ser.t <- FindClusters(ser.t,res=0.7)

DimPlot(ser.t,label=T)
FeaturePlot(ser.t,c("CD3D","CD8A","CD4","FCGR3A"))

## Add T cell and NK cell types
data.to.add <- vector("logical",length=ncol(ser.t))

data.to.add[ser.t@meta.data$RNA_snn_res.0.7=="0" |
              ser.t@meta.data$RNA_snn_res.0.7=="1" |
              ser.t@meta.data$RNA_snn_res.0.7=="2" |
              ser.t@meta.data$RNA_snn_res.0.7=="4" |
              ser.t@meta.data$RNA_snn_res.0.7=="5" |
              ser.t@meta.data$RNA_snn_res.0.7=="7" |
              ser.t@meta.data$RNA_snn_res.0.7=="10" |
              ser.t@meta.data$RNA_snn_res.0.7=="11"] <- "CD4 T cells"

data.to.add[ser.t@meta.data$RNA_snn_res.0.7=="3" |
              ser.t@meta.data$RNA_snn_res.0.7=="8" |
              ser.t@meta.data$RNA_snn_res.0.7=="9" |
              ser.t@meta.data$RNA_snn_res.0.7=="12"] <- "CD8 T cells"

data.to.add[ser.t@meta.data$RNA_snn_res.0.7=="6" |
              ser.t@meta.data$RNA_snn_res.0.7=="13"] <- "NK cells"

ser.t[["cell_types"]] <- data.to.add

## Check out T cell types
DimPlot(ser.t,group.by="cell_types")

## Add back to overall object
ser@meta.data$cell_types[match(colnames(ser.t),colnames(ser))] <- ser.t@meta.data$cell_types
DimPlot(ser,group.by="cell_types",split.by="sample_type")

## Metadata
overall.metadata <- ser@meta.data
format(object.size(overall.metadata),unit="MB")

## UMAP
overall.umap <- Embeddings(ser,reduction="umap")
format(object.size(overall.umap),unit="MB")

## Ligand/receptor matrix
load("data/ramilowski_pairs.rda")
ligs.recs.all <- unique(c(unique(as.character(ramilowski_pairs$ligand)),
                   unique(as.character(ramilowski_pairs$receptor))))
ligs.recs.use <- ligs.recs.all[ligs.recs.all %in% rownames(ser)]
overall.lig.rec <- GetAssayData(ser,slot="data",assay="RNA")[ligs.recs.use,]

# Filter to those expressed >1000 counts and <20000 counts
ligs.keep <- apply(overall.lig.rec,1,sum)[apply(overall.lig.rec,1,sum)>1000 &
                                            apply(overall.lig.rec,1,sum)<20000]
filtered.lig.rec <- overall.lig.rec[names(ligs.keep),]
format(object.size(filtered.lig.rec),unit="MB")

## Enforce consistent naming
overall_metadata <- overall.metadata
overall_umap <- overall.umap
filtered_lig_rec <- filtered.lig.rec

## Save minimum files
save(overall_metadata,file="data/overall_metadata.rda")
save(overall_umap,file="data/overall_umap.rda")
save(filtered_lig_rec,file="data/filtered_lig_rec.rda")
