#' Assess ligand and receptor interactions across groups of cells
#'
#' @param input_object Seurat object to create ligand receptor matrices from
#'
#' @param metadata_grouping Grouping variable in the Seurat object metadata used
#' to define the groups of cells. Default is "cell_types".
#'
#' @param ligand_receptor_pairs Data.frame of ligands, receptors and
#' interactions in the format of ramilowski_pairs provided by this package.
#' Defaults is "ramilowski_pairs".
#'
#' @param number_cells_required Number of cells per group required to perform
#' analysis of ligand/receptor interactions. Defaults to 100.
#'
#' @param min_expression Minimum expression in counts to consider a ligand or
#' receptor for interactions analysis. A sensible default is set to 1000, but is
#' dataset dependent. This is meant to filter out lowly expressed ligands and
#' receptors.
#'
#' @param max_expression Maxmium expression in counts to consider a ligand or
#' receptor for interactions analysis. A sensible default is set to 20000, but is
#' dataset dependent. This is meant to filter out ubiquitously expressed ligands
#' and receptors.
#'
#' @param scramble_times Number of times to scamble ligand/receptor interactions to
#' create a background distribution for statistical comparison.
#'
#' @return Comprehensive tibble of cognate ligand and receptor interactions and
#' statistical significance of these interactions across cell types.
#'
#' @export

celltalk <- function(input_object,
  metadata_grouping=cell_types,
  ligand_receptor_pairs=ramilowski_pairs,
  number_cells_required=100,
  min_expression=1000,
  max_expression=20000,
  scramble_times=10
) {

cell_types <- ramilowski_pairs <- NULL

  expr_mat_list <- create_expression_matrices(input_object,
    metadata_grouping,
    ligand_receptor_pairs,
    number_cells_required,
    min_expression,
    max_expression)

  lig_rec_means <- ligand_receptor_means(ligand_mean_dataframe=expr_mat_list[[1]],
  receptor_mean_dataframe=expr_mat_list[[2]],
  ligand_receptor_pairs)

  lig_rec_scram <- ligand_receptor_scramble(lig_rec_means,
  expr_mat_list[[1]],
  expr_mat_list[[2]],
  ligand_receptor_pairs,
  scramble_times)

  ligand_receptor_stats(lig_rec_means,lig_rec_scram)

}
