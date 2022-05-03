#' Statistically compare interactions between groups of samples
#'
#' @param seurat_object Seurat object containing expression data acorss all
#' samples and cells
#'
#' @param interaction_stats Tibble from celltalk function, usually filtered to
#' the top significant ligand and receptor interactions of interest
#'
#' @param sample_replicates Name of the meta.data column in a Seurat object that
#' has the samples of the individual replicate samples
#'
#' @param sample_group Name of the meta.data column in a Seurat object that
#' has the name of the sample group
#'
#' @param metadata_grouping Name of the meta.data column in a Seurat object that
#' has the name of the groups of cells to evaluate (e.g. "cell_types" containing
#' previously identified cell types)
#'
#' @return A list of linear models containing the coefficients and statistics for
#' assessing differences in joint ligand and receptor expression between sample
#' groups. Note that this requires replicate samples from each samples group.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @importFrom dplyr recode
#' @importFrom stats lm
#'
#' @export

compare_group_interactions <- function(seurat_object,
  interaction_stats,
  sample_replicates,
  sample_groups,
  metadata_grouping) {

  interactions_compare <- id_interactions(interaction_stats)

  extracted_sample_groups <- extract_sample_group_replicates(seurat_object,sample_groups)

  replicate_scores_frame <- cell_type_lig_rec_frame(interactions_compare,seurat_object,sample_replicates,sample_groups,metadata_grouping)
  # Replace dots with dashes in column names e.g. HLA.F should be HLA-F
  colnames(replicate_scores_frame) <- gsub("\\.","-",colnames(replicate_scores_frame))

  interaction_scores_frame <- cell_type_ligand_receptor_score(replicate_scores_frame,interactions_compare,
    extracted_sample_groups)

  input_formula <- as.formula(paste("scores~sample_groups",sep=""))

  mod_res <- interaction_scores_frame %>%
    mutate(unique_group=paste(cell_type1,ligand,cell_type2,receptor,sep="_")) %>%
    split(.$unique_group) %>%
    map(~lm(input_formula,data=.))

}
