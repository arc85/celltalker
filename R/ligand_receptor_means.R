#' Assess joint means between ligands and receptors across cell types
#'
#' @param ligand_mean_dataframe Mean expression values for ligands across cells
#' from a given condition
#'
#' @param receptor_mean_dataframe Mean expression values for receptors across
#' cells from a given condition
#'
#' @param ligand_receptor_pairs Data.frame of ligands, receptors and
#' interactions in the format of ramilowski_pairs provided by this package
#'
#' @param number_cell_types The number of unique cell types in a given sample
#' or condition. This will be used to assess ligand and receptor interactions
#' across cell types.
#'
#' @return Generates a tibble containing joint means for ligand receptor
#' interaction pairs across cell types.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#' @importFrom dplyr recode
#'
#' @export

ligand_receptor_means <- function(ligand_mean_dataframe,receptor_mean_dataframe,ligand_receptor_pairs,number_cell_types) {

# Bind variables
Var1 <- Var2 <- interaction_pairs <- NULL

lig.rec.pairs <- vector("list",length=nrow(ligand_receptor_pairs))

for (i in 1:nrow(ligand_receptor_pairs)) {

	lig.use <- as.character(ligand_receptor_pairs[i,"ligand"])
	rec.use <- as.character(ligand_receptor_pairs[i,"receptor"])

	lig.status <- !lig.use %in% rownames(ligand_mean_dataframe)
	rec.status <- !rec.use %in% rownames(receptor_mean_dataframe)

	if (any(lig.status,rec.status)) {

		lig.rec.pairs[[i]] <- NA

	} else {

		comparisons <- expand.grid(seq(1:number_cell_types),seq(1:number_cell_types))
		all.lig.rec <- vector("list")
		test.frame <- data.frame(ligand_mean_dataframe[lig.use,],receptor_mean_dataframe[rec.use,])

		for (q in 1:nrow(comparisons)) {

			all.lig.rec[[q]] <- mean(c(test.frame[comparisons[q,1],1],test.frame[comparisons[q,2],2]))

		}

		lig.rec.pairs[[i]] <- data.frame(do.call(rbind,all.lig.rec))
		colnames(lig.rec.pairs[[i]])[1] <- paste(lig.use,rec.use,sep="_")

	}

}

interact.frame <- do.call(cbind,lig.rec.pairs)
remove.cols <- apply(interact.frame,2,function(x) any(is.na(x)))
interact.frame <- interact.frame[,-which(remove.cols)]
remove.cols2 <- apply(interact.frame,2,function(x) all(x==0))
interact.frame <- interact.frame[,!remove.cols2]

level.key <- colnames(ligand_mean_dataframe)
names(level.key) <- seq(1,number_cell_types,1)

new.names <- comparisons %>%
	data.frame() %>%
	mutate(Var1=recode(Var1,!!!level.key)) %>%
	mutate(Var2=recode(Var2,!!!level.key))

rownames(interact.frame) <- paste(new.names[,1],new.names[,2],sep="_")

as_tibble(interact.frame,rownames="interaction_pairs")

}
