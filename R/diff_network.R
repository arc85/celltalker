#' Differential network analysis for putatitve ligand receptor interactions between groups
#'
#' Takes the output tibble from \code{putative_intearctions} and compares the networks of putative ligand receptor interactions between groups.
#'
#' @param putative.interactions.tib A tibble output from \code{putative_interactions}.
#' @param group1 Sample group to evaluate for unqiue ligands or receptors participating in interactions.
#' @param group2 Sample group as a comparitor to group1
#' @param interaction Name of a ligand_receptor interaction to evaluate
#'
#' @return Returns a tibble containing the differential graphs for group1 versus group2 and group2 versus group1, and the differential edges from group1 versus group2 and group2 versus group1. The differential edges are derived from the igraph package.
#'
#' @export


diff_network <- function(putative.interactions.tib,group1,group2,interaction) {

group1.res <- pull(putative.interactions.tib[putative.interactions.tib$name %in% group1,2])[[1]][[interaction]]
group2.res <- pull(putative.interactions.tib[putative.interactions.tib$name %in% group2,2])[[1]][[interaction]]

if (any(sapply(group2.res,function(x) is.null(x)))) {

	d <- expand.grid(group1.res)
	d[,1] <- as.character(d[,1])
	d[,2] <- as.character(d[,2])
	vertices <- unique(unlist(unname(group1.res)))
	net <- graph_from_data_frame(d=d,directed=T,vertices=vertices)

	enframe(list("diff_group1v2"=net))

} else if (any(sapply(group1.res,function(x) is.null(x)))) {

	d2 <- expand.grid(group2.res)
	d2[,1] <- as.character(d2[,1])
	d2[,2] <- as.character(d2[,2])
	vertices2 <- unique(unlist(unname(group2.res)))
	net2 <- graph_from_data_frame(d=d2,directed=T,vertices=vertices2)

	enframe(list("diff_group2v1"=net2))

} else {

d <- expand.grid(group1.res)
d[,1] <- as.character(d[,1])
d[,2] <- as.character(d[,2])
vertices <- unique(unlist(unname(group1.res)))
net <- graph_from_data_frame(d=d,directed=T,vertices=vertices)

d2 <- expand.grid(group2.res)
d2[,1] <- as.character(d2[,1])
d2[,2] <- as.character(d2[,2])
vertices2 <- unique(unlist(unname(group2.res)))
net2 <- graph_from_data_frame(d=d2,directed=T,vertices=vertices2)

g <- difference(net,net2)
rev.g <- difference(net2,net)

dif.g <- get.edgelist(g)
colnames(dif.g) <- c("group1","group2")
dif.rev.g <- get.edgelist(rev.g)
colnames(dif.rev.g) <- c("group1","group2")

enframe(list("diff_group1v2"=g,"diff_group2v1"=rev.g,"diff_edge_group1v2"=dif.g,"diff_edge_group2v1"=dif.rev.g))

}

}
