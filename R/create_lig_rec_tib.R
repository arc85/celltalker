#' Create a tibble of consistently expressed ligands and receptors
#'
#' @param exp.tib The nested tibble that was the output from \code{reshape_matrices}.
#' @param clusters A named vector containing the cluster identities of each cell as a factor, where the names match the column names from the count matrix.
#' @param groups A named vector containing the group identities of each cell as a factor, where the names match the columna names from the count matrix.
#' @param replicates A named vector containing the replicate identities of each sample as a factor, where the names match the column names from the count matrix. Assumes there are biological replicates of each sample.
#' @param cells.reqd An integer value representing the number of cells required to express a given ligand or receptor within a cluster from an individual replicate.
#' @param freq.pos.reqd A fraction from 0 to 1 representing the fraction of individual replicate samples that must express a ligand or receptor above the \code{cells.reqd} threshold to be counted as a ligand or receptor in a cluster from a given group.
#' @param ligands.and.receptors A matrix with 3 columns, the first of which corresponds to ligands, the second of which corresponds to receptors and the last of which corresponds to ligand_receptor interactions. Column names should be "ligand", "receptor", and "pairs" respectively.
#'
#' @return Returns a nested tibble containing the ligands and receptors consistently expressed within clusters from each group.
#'
#' @export


create_lig_rec_tib <- function(exp.tib,clusters,groups,replicates,cells.reqd,freq.pos.reqd,ligands.and.receptors) {

clusters <- as.factor(clusters)

replicate.tab <- vector("list",length=length(levels(as.factor(groups))))
names(replicate.tab) <- levels(as.factor(groups))

for (a in 1:length(levels(as.factor(groups)))) {

pid.layer <- unnest(exp.tib[a,2],cols="samples")
pid.cell.num <- vector("list",length=nrow(pid.layer))

lig.rec.res <- vector("list",length=length(levels(clusters)))
names(lig.rec.res) <- levels(clusters)

for (z in 1:length(levels(clusters))) {

	for (i in 1:nrow(pid.layer)) {
		pid.cell.num[[i]] <- Reduce(rbind,unnest(pid.layer[i,2],cols="expr.matrices") %>% transmute(n.rows=map(expr.matrices,nrow)) %>% pull(n.rows))
	}

	n.cells.cluster <- Reduce(rbind,lapply(pid.cell.num,function(x) x[z,]))

	if (all(n.cells.cluster<cells.reqd)) {

		lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)

	} else if (sum(n.cells.cluster>cells.reqd)/nrow(pid.layer)<=freq.pos.reqd) {

		lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)

		} else {

		cluster.pos <- pid.layer %>% transmute(sample,pos=map(expr.matrices,~.x[[z]])) %>% pull(pos)
		cluster.pos <- cluster.pos[n.cells.cluster>cells.reqd]
		cols.to.drop <- c("cell.names","replicate.id","group.id","cluster.id")
		cluster.pos <- lapply(cluster.pos,function(x) select(x,-cols.to.drop))

		genes.pos <- lapply(cluster.pos,function(x) apply(x,2,function(y) sum(y>0)>freq.pos.reqd))
		genes.pos.vec <- Reduce(rbind,genes.pos)

		if (is.null(dim(genes.pos.vec))) {
			if (1/nrow(pid.layer)>freq.pos.reqd) {
			genes.to.include <- names(genes.pos.vec)[genes.pos.vec]
		} else {lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)}
		} else {

		genes.to.include <- apply(genes.pos.vec,2,function(x) sum(x)/nrow(pid.layer)>freq.pos.reqd)
		genes.to.include <- names(genes.to.include)[genes.to.include]

		ligands <- genes.to.include[genes.to.include %in% ligands.and.receptors$ligand]
		receptors <- genes.to.include[genes.to.include %in% ligands.and.receptors$receptor]

		lig.rec.res[[z]] <- list(ligands=ligands,receptors=receptors)

		}

	}

}

replicate.tab[[a]] <- tibble(cluster.id=levels(clusters),ligands.and.receptors=lig.rec.res)

}

lig.rec.tab <- enframe(replicate.tab,name="group",value="lig.rec.exp")

}
