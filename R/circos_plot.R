#' Creates a circos plot from the list of ligands and receptors
#'
#' @param interactions1 List of ligands and receptors and their associated cell types (i.e. the output from one sample group from the putative_interactions function). If interactions2 is NULL, then this will only plot interactions found in the list of ligands and receptors from the sample group specified.
#'
#' @param interactions2 List of ligands and receptors and their associated cell types (i.e. the output from a second sample group from the putative_interactions function). Default is NULL, but if this value is non-null, will plot both sets of interactions present in "interactions1" and "interactions2". This group is considered the "control" group if present, and unique ligand/receptor interactions from "interactions1" will be plotted.
#'
#' @param ligand.col Color of the bars defining the ligands
#'
#' @param receptor.col Color of the bars definig the receptors
#'
#' @param interactions.col Color of the links connecting ligands and receptors. Default is light gray
#'
#' @param unique.vs.control.interactions.col Optional color of the links connecting the ligand and receptor interactions that are unique between the sample group and the control
#'
#' @return Generates a circos plot connecting ligands and receptors across cell types for a given sample group
#'
#' @export

circos_plot <- function(interactions1,interactions2=NULL,clusters,ligand.col="#C8A757",receptor.col="#70BC71",interactions.col="lightgray",unique.vs.control.interactions.col="darkorange") {

  clusters <- levels(clusters)

  if (is.null(interactions2)) {

  interactions.list <- list(interactions1)

  } else {

    interactions.list <- list(interactions1,interactions2)

  }

	#Setup overall circos plot

  if (length(interactions.list)==1) {

  non.neg.int <- sapply(sapply(interactions.list,function(y) y),function(x) !is.null(x[[1]]))
	shared.int <- interactions.list[[1]][non.neg.int]

  ligands <- unique(sapply(strsplit(names(shared.int),split="_"),`[[`,1))
  receptors <- unique(sapply(strsplit(names(shared.int),split="_"),`[[`,2))

  } else {

    all.int <- c(names(interactions.list[[1]]),names(interactions.list[[2]]))
    non.neg.int <- sapply(sapply(interactions.list,function(y) y),function(x) !is.null(x[[1]]))
    all.int <- unique(all.int[non.neg.int])

    group1.pos <- sapply(interactions.list[[1]],function(x) !is.null(x[[1]]))
    group2.pos <- sapply(interactions.list[[2]],function(x) !is.null(x[[1]]))

    shared.names <- intersect(names(interactions.list[[1]])[group1.pos],names(interactions.list[[2]])[group2.pos])

    group1.shared <- interactions.list[[1]][shared.names]
    group2.shared <- interactions.list[[2]][shared.names]
    shared.int1 <- mapply(union,sapply(group1.shared,function(x) x[[1]]),sapply(group2.shared,function(x) x[[1]]),SIMPLIFY=FALSE)
    shared.int2 <- mapply(union,sapply(group1.shared,function(x) x[[2]]),sapply(group2.shared,function(x) x[[2]]),SIMPLIFY=FALSE)

    shared.int <- vector("list",length=length(shared.names))
    for (i in 1:length(shared.int)) {
      shared.int[[i]] <- list("ligand.cells"=shared.int1[[i]],"receptor.cells"=shared.int2[[i]])
    }
    names(shared.int) <- shared.names

    group1.uni.names <- setdiff(names(interactions.list[[1]])[group1.pos],names(interactions.list[[2]])[group2.pos])
    group1.uni <- interactions.list[[1]][group1.uni.names]

    ligands <- unique(sapply(strsplit(all.int,split="_"),`[[`,1))

    receptors <- unique(sapply(strsplit(all.int,split="_"),`[[`,2))

  }

	circos.data <- data.frame(matrix(data=0,nrow=length(clusters)*length(c(ligands,receptors)),ncol=3))
	colnames(circos.data) <- c("cell.type","genes","ligand.rec")

	num.lig.rec <- sum(length(ligands),length(receptors))

	seq1 <- seq(1,num.lig.rec*length(clusters),num.lig.rec)
	seq2 <- seq(num.lig.rec,num.lig.rec*length(clusters)+1,num.lig.rec)

	for (i in 1:length(seq2)) {
		circos.data[seq1[i]:seq2[i],1] <- clusters[i]
		circos.data[seq1[i]:seq2[i],2] <- c(ligands,receptors)
		circos.data[seq1[i]:seq2[i],3] <- c(rep("lig",length(ligands)),rep("rec",length(receptors)))
	}

	circos.data$ref.num <- seq(1,num.lig.rec*length(clusters),1)


	#Plot circos

	circos.par("cell.padding"=c(0.02,0,0.02,0),"track.height"=0.1)
	circos.initialize(factors=circos.data$cell.type,x=seq(1,nrow(circos.data),1))

	circos.track(factors=circos.data$clusters,ylim=c(0,1),
		panel.fun=function(x,y) {
			circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2] + uy(5,"mm"),CELL_META$sector.index)
			circos.rect(xleft=CELL_META$xlim[1],ybottom=CELL_META$cell.ylim[1],xright=CELL_META$xlim[1]+length(ligands),ytop=CELL_META$cell.ylim[2],col=ligand.col)
			circos.rect(xleft=CELL_META$xlim[1]+length(ligands),ybottom=CELL_META$cell.ylim[1],xright=CELL_META$xlim[2],ytop=CELL_META$cell.ylim[2],col=receptor.col)
		}
	)


	#Draw links - shared

	for (i in 1:length(names(shared.int))) {
		link <- names(shared.int)[i]

		lig.link <- sapply(strsplit(link,split="_"),`[[`,1)
		rec.link <- sapply(strsplit(link,split="_"),`[[`,2)

		link.1 <- shared.int[[i]][["ligand.cells"]]
		link.2 <- shared.int[[i]][["receptor.cells"]]

		for (a in 1:length(link.1)) {

		circos.sub1 <- subset(circos.data,cell.type==link.1[a])

			for (z in 1:length(link.2)) {

				circos.sub2 <- subset(circos.data,cell.type==link.2[z])

				ref.num1 <- circos.sub1$ref.num[which(circos.sub1$genes==lig.link)]
				ref.num2 <- circos.sub2$ref.num[which(circos.sub2$genes==rec.link)]

				circos.link(link.1[a],ref.num1,link.2[z],ref.num2,col=interactions.col,h.ratio=0.2)
		}

	}
}

if (length(interactions.list)==1) {
  print("Plotted interactions")
} else {
  print("Plotted shared interactions")
}

if (length(interactions.list)>1) {

#Draw links - unique

for (i in 1:length(names(group1.uni))) {
  link <- names(group1.uni)[i]

  lig.link <- sapply(strsplit(link,split="_"),`[[`,1)
  rec.link <- sapply(strsplit(link,split="_"),`[[`,2)

  link.1 <- group1.uni[[i]][["ligand.cells"]]
  link.2 <- group1.uni[[i]][["receptor.cells"]]

  for (a in 1:length(link.1)) {

  circos.sub1 <- subset(circos.data,cell.type==link.1[a])

    for (z in 1:length(link.2)) {

      circos.sub2 <- subset(circos.data,cell.type==link.2[z])

      ref.num1 <- circos.sub1$ref.num[which(circos.sub1$gene==lig.link)]
      ref.num2 <- circos.sub2$ref.num[which(circos.sub2$gene==rec.link)]

      circos.link(link.1[a],ref.num1,link.2[z],ref.num2,col=unique.vs.control.interactions.col,h.ratio=0.2)
  }

}
}

print("Plotted unique interactions")

}

}
