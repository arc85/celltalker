#' Evaluate differential interactions
#'
#' Evaluates whether a given cluster participates in putative ligand/receptor interactions between groups. Can search for a cluster with a ligand or receptor participating in differential interactions.

#' @param putative.interactions.tib A tibble output from \code{putative_interactions}.
#' @param group1 Sample group to evaluate for unqiue ligands or receptors participating in interactions.
#' @param group2 Sample group as a comparitor to group1
#' @param cluster A cluster to evaluate unique ligand/receptor interactions between groups.
#' @param search.ligand Logical dictating whether to search for the presence of a unique ligand between groups. Must be TRUE if \code{search.receptor} is FALSE.
#' @param searh.receptor Logical dictating whether to search for the presence of a unique ligand between groups. Must be TRUE if \code{search.ligand} is FALSE.
#' @param ligands.and.receptors A matrix with 3 columns, the first of which corresponds to ligands, the second of which corresponds to receptors and the last of which corresponds to ligand_receptor interactions. Column names should be "ligand", "receptor", and "pairs" respectively.
#'
#' @return Returns a tibble containing the interactions where a cluster has a unique ligand or receptor interaction.
#'
#' @export

diff_interactions <- function(putative.interactions.tib,group1,group2,cluster,search.ligand=TRUE,search.receptor=FALSE,ligands.and.receptors) {

	if (search.ligand==TRUE & search.receptor==TRUE) {
		return(print("Error: Only one of 'search.ligand' or 'search.receptor' can be TRUE"))
	}
	else if (search.ligand==TRUE) {
		search.index <- 1
	} else if (search.receptor==TRUE) {
		search.index <- 2
	} else {
		return(print("One of 'search.ligand' or 'search.receptor' must be TRUE"))
	}

	set1 <- unnest(putative.interactions.tib[putative.interactions.tib$name %in% group1,2]) %>% pull()
	set2 <- unnest(putative.interactions.tib[putative.interactions.tib$name %in% group2,2]) %>% pull()

	set1.search <- lapply(set1,function(x) x[search.index])
	set2.search <- lapply(set2,function(x) x[search.index])

	set1.cluster <- grep(cluster,set1.search)
	set2.cluster <- grep(cluster,set2.search)

	set1.pairs <- ligands.and.receptors$pair[set1.cluster]
	set2.pairs <- ligands.and.receptors$pair[set2.cluster]

	set1.diff <- setdiff(set1.pairs,set2.pairs)
	set2.diff <- setdiff(set2.pairs,set1.pairs)

	set.diff.list <- list(set1.diff,set2.diff)
	names(set.diff.list) <- c(group1,group2)

	enframe(set.diff.list,name="group",value="unique_interactions")

}
