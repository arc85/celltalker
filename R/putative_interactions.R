#' Create a tibble of consistently expressed ligands and receptors
#'
#' Takes the tibble output from \code{create_lig_rec_tib} and compares ligand and receptor expression withn and between clusters within a given sample group. Uses interactions described in \code{ligands.and.receptors}.
#'
#' @param ligand.receptor.tibble The nested tibble that was the output from \code{create_lig_rec_tib}.
#' @param clusters A named vector containing the cluster identities of each cell as a factor, where the names match the column names from the count matrix.
#' @param groups A named vector containing the group identities of each cell as a factor, where the names match the columna names from the count matrix.
#' @param freq.group.in.cluster A fraction from 0 to 1 describing the frequency of cells from a given group required to be in a cluster to consider that cluster for interactions.
#' @param ligands.and.receptors A matrix with 3 columns, the first of which corresponds to ligands, the second of which corresponds to receptors and the last of which corresponds to ligand_receptor interactions. Column names should be "ligand", "receptor", and "pair" respectively.
#'
#' @return Returns a nested tibble containing putative interactions between clusters for a given group, where each putative interaction is a member of a list and the clusters expressing cognate ligands and receptore are in the list.
#'
#' @export


putative_interactions <- function(ligand.receptor.tibble,clusters,groups,freq.group.in.cluster,ligands.and.receptors) {

group.list <- vector("list",length=length(levels(as.factor(groups))))
names(group.list) <- levels(as.factor(groups))

clusters.to.include <- vector("list",length=length(group.list))
names(clusters.to.include) <- names(group.list)
clusters.per.group <- table(clusters,groups)/rowSums(table(clusters,groups))

for (q in 1:length(group.list)) {

	interactions.list <- vector("list",length=nrow(ligands.and.receptors))
	names(interactions.list) <- ligands.and.receptors$pair

	sub.list <- list("ligand.cells"=NULL,"receptor.cells"=NULL)

	for (i in 1:length(interactions.list)) {
		interactions.list[[i]] <- sub.list
		}

group.unnest <- unnest(ligand.receptor.tibble[q,2])
clusters.to.use <- rownames(clusters.per.group)[clusters.per.group[,q]>freq.group.in.cluster]
group.unnest <- group.unnest[group.unnest$cluster.id %in% clusters.to.use,]

for (z in 1:nrow(group.unnest)) {

cell.ligs <- pull(group.unnest[z,2])[[1]]$ligands

rec.list <- lapply(pull(group.unnest),function(x) x[[2]])
names(rec.list) <- group.unnest$cluster.id

for (a in 1:length(rec.list)) {

	if (is.null(rec.list[[a]])) {

	} else {

expanded <- expand.grid(cell.ligs,rec.list[[a]])
expanded$pair <- paste(expanded[,1],expanded[,2],sep="_")
interactions <- expanded$pair[expanded$pair %in% ramilowski_pairs$pair]

}

if (length(interactions)==0) {

} else {

for (i in 1:length(interactions.list[interactions])) {
	interactions.list[interactions][[i]]$ligand.cells <- c(interactions.list[interactions][[i]]$ligand.cells,names(rec.list)[z])
	interactions.list[interactions][[i]]$ligand.cells <- unique(interactions.list[interactions][[i]]$ligand.cells)
	interactions.list[interactions][[i]]$receptor.cells <- c(interactions.list[interactions][[i]]$receptor.cells,names(rec.list)[a])
	interactions.list[interactions][[i]]$receptor.cells <- unique(interactions.list[interactions][[i]]$receptor.cells)
}

}

}

}

group.list[[q]] <- interactions.list

}

group.tab <- enframe(group.list)

}
