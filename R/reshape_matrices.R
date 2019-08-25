#' Reshape matrices
#'
#' Takes an overall count matrix and metadata describing the clusters, groups, replicates and ligands and receptors as input and generates individual gene expression matrices for each cluster within each replicate within each group. Returns a nested tibble with the individual expression matrices.
#'
#' @param count.matrix A raw count matrix with rows as genes and columns as cells.
#' @param clusters A named vector containing the cluster identities of each cell as a factor, where the names match the column names from the count matrix.
#' @param groups A named vector containing the group identities of each cell as a factor, where the names match the columna names from the count matrix.
#' @param replicates A named vector containing the replicate identities of each sample as a factor, where the names match the column names from the count matrix. Assumes there are biological replicates of each sample.
#' @param ligands.and.receptors A matrix with 3 columns, the first of which corresponds to ligands, the second of which corresponds to receptors and the last of which corresponds to ligand_receptor interactions. Column names should be "ligand", "receptor", and "pair" respectively.
#'
#' @return Returns a nested tibble containing the individual expression matrices, with the first level as the groups and the second level as the individual replciates within each group.
#'
#' @export

reshape_matrices <- function(count.matrix,clusters,groups,replicates,ligands.and.receptors) {

	#Filter by ligands and receptors
	union.lig.rec <- union(ligands.and.receptors$ligand,ligands.and.receptors$receptor)
	mat.fil <- count.matrix[rownames(count.matrix) %in% union.lig.rec,]

	#Create combined metadata
	comb.meta <- data.frame("replicate.id"=replicates,"group.id"=groups,"cluster.id"=clusters)

	#Add metadata to count matrix
	mat.tib <- data.frame(as.matrix(mat.fil))
	mat.tib <- t(mat.tib)
	mat.tib <- cbind(mat.tib,comb.meta)
	mat.tib <- rownames_to_column(mat.tib,var="cell.names")
	mat.tib <- as_tibble(mat.tib)

	#Create nested data.frames by splitting and combining
	sp.rep.id <- mat.tib %>% split(.$replicate.id)
	sp.clust.id <- sp.rep.id %>% map(~.x %>% split(.$cluster.id))
	sp.clust.id.nest <- enframe(sp.clust.id,'sample')

	ref.table <- table(groups,replicates)

	fun <- function(x) { names(which.max(ref.table[,colnames(ref.table) %in% sp.clust.id.nest$sample[x]])) }

	res <- data.frame("replicate.id"=sp.clust.id.nest$sample,"group"=NA)
	for (i in 1:ncol(ref.table)) {
		res[i,2] <- fun(i)
	}

	groups.list <- res %>% split(.$group)

	listing <- vector('list',length=length(groups.list))
	for (i in 1:length(groups.list)) {
		listing[[i]] <- sp.clust.id.nest[sp.clust.id.nest$sample %in% groups.list[[i]]$replicate.id,]
	}
	names(listing) <- names(groups.list)

	groups.samples.nested <- enframe(listing)

}
