#' Create a boxplot of the joint mean score for a specific ligand receptor
#' interaction between two groups of samples
#'
#' @param ligand Ligand of interest
#'
#' @param receptor Receptor of interest
#'
#' @param cell_type1 Cell type expressing the ligand of interest
#'
#' @param cell_type2 Cell type expressiong the receptor of interest
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
#' @return A ggplot boxplot of the joint mean for a ligand/receptor interaction
#' in each sample group.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter theme_bw xlab ylab guides guide_legend
#'
#' @export


boxplot_group_interaction <- function(ligand,receptor,cell_type1,cell_type2,sample_replicates,sample_group_1,sample_group_2) {

  score1 <- cell_type_score(ligand=ligand,
    receptor=receptor,
    cell_type1=cell_type1,
    cell_type2=cell_type2,
    sample_replicates=sample_replicates,
    sample_group=sample_group_1
  )

  score2 <- cell_type_score(ligand=ligand,
    receptor=receptor,
    cell_type1=cell_type1,
    cell_type2=cell_type2,
    sample_replicates=sample_replicates,
    sample_group=sample_group_2
  )

  extract_names1 <- sample_group_1[["sample_type"]][1,1]
  extract_names2 <- sample_group_2[["sample_type"]][1,1]

  res <- data.frame(sample_id=c(names(score1),names(score2)),
  sample_group=c(rep(extract_names1,length=length(score1)),
    rep(extract_names2,length=length(score2))),
    Interaction_scores=c(score1,score2))

  ggplot(res,aes(x=sample_group,y=Interaction_scores,colour=sample_group)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1,height=0) +
    theme_bw() +
    xlab("Sample group") +
    ylab("Interaction scores") +
    guides(colour=guide_legend(title="Sample group"))

}
