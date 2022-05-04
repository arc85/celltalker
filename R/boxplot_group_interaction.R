#' Create a boxplot of the joint mean score for a specific ligand receptor
#' interaction between two groups of samples
#'
#' @param seurat_object Seurat object containing expression levels across all
#' cells and sample types
#'
#' @param interaction_stats Tibble of interaction statistics from a sample group
#'
#' @param interaction_stats Tibble from celltalk function, usually filtered to
#' the top significant ligand and receptor interactions of interest
#'
#' @param sample_replicates Name of the meta.data column in a Seurat object that
#' has the samples of the individual replicate samples
#'
#' @param sample_groups Name of the meta.data column in a Seurat object that
#' has the name of the sample group
#'
#' @param metadata_grouping Name of the meta.data column in a Seurat object that
#' has the name of the groups of cells to evaluate (e.g. "cell_types" containing
#' previously identified cell types)
#'
#' @param ligand Name of the ligand in the ligand/receptor pair of interest
#'
#' @param receptor Name of the receptor in the ligand/receptor pair of interest
#'
#' @param cell_type1 Name of the cell type expressing the ligand in the ligand/
#' receptor pair of interest
#'
#' @param cell_type2 Name of the cell type expressing the receptor in the ligand/
#' receptor pair of interest
#'
#' @return A ggplot boxplot of the joint mean for a ligand/receptor interaction
#' in each sample group.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter theme_bw ylab
#'
#' @export

boxplot_group_interaction <- function(seurat_object,
  interaction_stats,
  sample_replicates,
  sample_groups,
  metadata_grouping,
  ligand,
  receptor,
  cell_type1,
  cell_type2) {

    `Sample group` <- scores <- `.` <- NULL

    interactions_compare <- id_interactions(interaction_stats)

    extracted_sample_groups <- extract_sample_group_replicates(seurat_object,sample_groups)

    replicate_scores_frame <- cell_type_lig_rec_frame(interactions_compare,seurat_object,sample_replicates,sample_groups,metadata_grouping)
    # Replace dots with dashes in column names
    colnames(replicate_scores_frame) <- gsub("\\.","-",colnames(replicate_scores_frame))

    interaction_scores_frame <- cell_type_ligand_receptor_score(replicate_scores_frame,interactions_compare,
      extracted_sample_groups)

    interaction_scores_frame %>%
      filter(cell_type1=={{cell_type1}}) %>%
      filter(cell_type2=={{cell_type2}}) %>%
      filter(ligand=={{ligand}}) %>%
      filter(receptor=={{receptor}}) %>%
      mutate(`Sample group`=sample_groups) %>%
      ggplot(.,aes(x=`Sample group`,y=scores,colour=`Sample group`)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(width=0.1,height=0) +
        theme_bw() +
        ylab("Joint ligand/receptor interaction")

}
