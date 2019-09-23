#' Creates a circos plot from the list of ligands and receptors
#'
#' @param interactions.list The output from the interactions_list_wrapper. Consists of a list of ligands and receptors in each cell type from each sample group
#'
#' @param sample.name Name of the sample group to plot from the interactions.list
#'
#' @param control.group Optional name of the control group to compare interactions in the sample group to. Default is false
#'
#' @param ligand.col Color of the bars defining the ligands
#'
#' @param receptor.col Color of the bars definig the receptors
#'
#' @param interactions.col Color of the links connecting ligands and receptors. Default is light gray
#'
#' @param unique.vs.control.interactions.col Optional color of the links connecting the ligand and receptor interactions that are unique between the sample group and the control
#'
#' @param metadata.matrix Dataframe containing the metadata from the expression dataset
#'
#' @return Generates a circos plot connecting ligands and receptors across cell types for a given sample group
#'
#' @export

circos_plot <- function(interactions1,interactions2=NULL,clusters,ligand.col="#C8A757",receptor.col="#70BC71",interactions.col="lightgray",unique.vs.control.interactions.col="darkorange") {

  if (is.null(interactions2)) {

  interactions.list <- interactions1

  } else {

    interactions.list <- list(interactions1,interactions2)

  }

	#Setup overall circos plot

  non.neg.int <- sapply(interactions.list,function(x) !is.null(x[[1]]))
	all.lig.rec.interactions <- interactions.list[non.neg.int]

	ligands <- unique(sapply(strsplit(names(all.lig.rec.interactions),split="_"),`[[`,1))

	receptors <- unique(sapply(strsplit(names(all.lig.rec.interactions),split="_"),`[[`,2))

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

	#Draw links

	for (i in 1:length(names(all.lig.rec.interactions))) {
		link <- names(all.lig.rec.interactions)[i]

		lig.link <- sapply(strsplit(link,split="_"),`[[`,1)
		rec.link <- sapply(strsplit(link,split="_"),`[[`,2)

		link.1 <- all.lig.rec.interactions[[i]][["ligand.cells"]]
		link.2 <- all.lig.rec.interactions[[i]][["receptor.cells"]]

		for (a in 1:length(link.1)) {

		circos.sub1 <- subset(circos.data,cell.type==link.1[a])

			for (z in 1:length(link.2)) {

				circos.sub2 <- subset(circos.data,cell.type==link.2[z])

				ref.num1 <- circos.sub1$ref.num[which(circos.sub1$gene==lig.link)]
				ref.num2 <- circos.sub2$ref.num[which(circos.sub2$gene==rec.link)]

				circos.link(link.1[a],ref.num1,link.2[z],ref.num2,col=interactions.col,h.ratio=0.2)
		}

	}
}

print("Plotted interactions")

}
