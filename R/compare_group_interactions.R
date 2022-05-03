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

compare_group_interactions <- function(seurat_object,
  interaction_stats,
  sample_replicates,
  sample_groups,
  p_value_filter=NULL,
  fdr_value_filter=0.05) {

  input_formula <- as.formula(paste("scores~sample_groups",sep=""))
  interactions_compare <- id_interactions(interaction_stats)

  extracted_sample_groups <- extract_sample_group_replicates()

  replicate_scores_frame <- cell_type_lig_rec_frame(interactions_compare,seurat_object,sample_replicates,sample_groups,cell_types)
  # Replace dots with dashes in column names
  colnames(replicate_scores_frame) <- gsub("\\.","-",colnames(replicate_scores_frame))

  ## Need to derive a way to extract sample groups
  interaction_scores_frame <- cell_type_ligand_receptor_score(replicate_scores_frame,interactions_compare,
    extracted_sample_groups)

  mod_res <- interaction_scores_frame %>%
    mutate(unique_group=paste(cell_type1,ligand,cell_type2,receptor,sep="_")) %>%
    split(.$unique_group) %>%
    map(~lm(input_formula,data=.))

  p_vals <- lapply(mod_res,function(x) {
    if (class(x) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(x)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  })

  fdr_vals <- p.adjust(p_vals,method="fdr")

  ## Filter results
  lapply(mod_res[fdr_vals<0.05],function(x) summary(x))

}
