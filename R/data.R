#' Expressed ligands and receptors.
#'
#' A dataset containing the filtered expressed of ligands and receptors from
#' a series of peripheral blood from 5 healthy donors and tonsil from 5
#' healthy donors undergoing tonsilectomy.
#'
#' @format A sparse matrix with 192 genes (rows) and 23734 cells (columns)
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324/}
"filtered_lig_rec"


#' UMAP of cells from 5 health blood donors and 5 healthy tonsil donors.
#'
#' The UMAP_1 and UMAP_2 coordinates of cells resulting from dimensionality
#' reduction of all cells from 5 healthy peripheral blood donors and 5 tonsil
#' tissue donors.
#'
#' @format A matrix with 23734 cells (rows) and 2 dimensions of UMAP embedding
#' (columns)
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324/}
"overall_umap"


#' Metadata from 5 health blood donors and 5 healthy tonsil donors.
#'
#' Overall metadata from analysis of 5 healthy blood donors and 5 healthy
#' tonsil tissue donors.
#'
#' @format A data.frame with 23734 cells (rows) and 8 variables:
#' \describe{
#'   \item{orig.ident}{name assigned by default to the Seurat object}
#'   \item{nCount_RNA}{number of RNA molecules counted per cell (1038--74386)}
#'   \item{nFeature_RNA}{number of RNA molecules counted per cell (41--6582)}
#'   \item{sample_id}{individual donor samples, HD_PBMC 1 to 5 and HD_Tonsil 1 to 5}
#'   \item{sample_type}{tissue of origin, either PBMC or Tonsil}
#'   \item{RNA_snn_res.0.3}{Louvian clustering results, clusters 0 to 13}
#'   \item{seurat_clusters}{clusters assigned by default, same as RNA_snn_res.0.3}
#'   \item{cell_types}{cell types identified by canonical gene expression (B cells,
#'    CD14-CD16+ monocytes, CD14+CD16- monocytes, CD1C+ DCs, CD4 T cells,
#'    CD8 T cells, NK cells, pDCs, Plasmablasts, RBCs)}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324/}
"overall_metadata"


#' List of known ligands, receptors and their iteractions.
#'
#' This is a curated list of ligands, receptors and their potential interactions.
#' The list is derived from a manuscript from Ramilowski et al. There are other
#' potential lists that could be used in place of or in addition to this list.
#' Notably, the CellPhoneDB (https://github.com/ventolab/CellphoneDB) package
#' has an extensively curated list including multiple receptor components.
#'
#' @format A data.frame with 2557 rows and 3 variables:
#' \describe{
#'   \item{ligand}{expressed ligands, total of 708 unique ligands}
#'   \item{receptor}{expressed receptors, total of 691 receptors}
#'   \item{pair}{interaction between and ligand and receptor, total of 2557
#'    unique interactions}
#' }
#' @source \url{https://www.nature.com/articles/ncomms8866}
"ramilowski_pairs"
