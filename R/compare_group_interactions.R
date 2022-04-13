#' Statistically compare interactions between groups of samples
#'
#' @param interaction_stats Tibble from celltalk function, usually filtered to
#' the top significant ligand and receptor interactions of interest
#'
#' @param sample_replicates Name of the meta.data column in a Seurat object that
#' has the samples of the individual replicate samples
#'
#' @param sample_group_1 The group of samples from which interactions will be
#' identified for statistical comparison
#'
#' @param sample_group_2 Comparitor group from which interactions from
#' sample_group_1 will be measured against
#'
#' @return A data.frame contain ligands, receptors, cell types and the
#' statistical significant of the joint means in sample_group_1 versus
#' the joint means in sample_group_2
#'
#' @export

compare_group_interactions <- function(interactions_stats,sample_replicates,sample_group_1,sample_group_2) {

  interactions_compare <- id_interactions(interactions_stats)

  stats_vector <- vector("list",length=nrow(interactions_compare))

  for (i in 1:nrow(interactions_compare)) {

    score1 <- cell_type_score(ligand=interactions_compare[i,"ligand"],
      receptor=interactions_compare[i,"receptor"],
      cell_type1=interactions_compare[i,"cell_type1"],
      cell_type2=interactions_compare[i,"cell_type2"],
      sample_replicates=sample_replicates,
      sample_group=sample_group_1
    )

    score2 <- cell_type_score(ligand=interactions_compare[i,"ligand"],
      receptor=interactions_compare[i,"receptor"],
      cell_type1=interactions_compare[i,"cell_type1"],
      cell_type2=interactions_compare[i,"cell_type2"],
      sample_replicates=sample_replicates,
      sample_group=sample_group_2
    )

    stats_vector[[i]] <- data.frame(mean_group1=mean(score1),
      mean_group2=mean(score2),
      p_value=wilcox.test(score1,score2)$p.value)

  }

  stats_frame <- do.call(rbind,stats_vector)

  cbind(interactions_compare,stats_frame)

}
