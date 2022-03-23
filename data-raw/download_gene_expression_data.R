## Download samples from Gene Expression Omnibus
## GSE139324
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324
## Dec 22 2021

## Identify samples from GEO
gds <- GEOquery::getGEO("GSE139324")
samples.use <- c("HD_PBMC_1","HD_PBMC_2","HD_PBMC_3","HD_PBMC_4","HD_PBMC_5",
  "HD_Tonsil_1","HD_Tonsil_2","HD_Tonsil_3","HD_Tonsil_4","HD_Tonsil_5")
data.use <- Biobase::phenoData(gds[[1]])@data[Biobase::phenoData(gds[[1]])@data$title %in% samples.use,]

## Create overall directory, subdirectories, and download from GEO
for (i in 1:nrow(data.use)) {

  dir.name <- paste("data-raw",data.use$title[i],sep="/")
  dir.create(dir.name)

  tmp.dwnload1 <- as.character(data.use[i,which(grepl("supplementary_file",colnames(data.use)))[1]])
  tmp.dwnload2 <- as.character(data.use[i,which(grepl("supplementary_file",colnames(data.use)))[2]])
  tmp.dwnload3 <- as.character(data.use[i,which(grepl("supplementary_file",colnames(data.use)))[3]])

  download.file(tmp.dwnload1,destfile=paste(dir.name,"barcodes.tsv.gz",sep="/"))
  download.file(tmp.dwnload2,destfile=paste(dir.name,"features.tsv.gz",sep="/"))
  download.file(tmp.dwnload3,destfile=paste(dir.name,"matrix.mtx.gz",sep="/"))

}
