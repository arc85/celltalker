#' Create ligand receptor matrices from a Seurat Object
#'
#' @param input_object Seurat object to create ligand receptor matrices from
#'
#' @param metadata_grouping Grouping variable in the Seurat object metadata used
#' to define the groups of cells. Default is "cell_types".
#'
#' @param ligand_receptor_pairs Data.frame of ligands, receptors and
#' interactions in the format of ramilowski_pairs provided by this package.
#' Defaults is "ramilowski_pairs".
#'
#' @return Generates a list, with the first element being the mean expression of
#' ligands in rows and cell groups in columns  and the second element being the
#' mean expression of receptors in rows and cell groups in columns.
#'
#' @keywords internal
#' @noRd

create_expression_matrices <- function(input_object,
  metadata_grouping,
  ligand_receptor_pairs,
  number_cells_required,
  min_expression,
  max_expression
  ) {

  # Define cells to keep
  cell_types_keep <- names(table(input_object[[metadata_grouping]]))[table(input_object[[metadata_grouping]])>number_cells_required]

  # Filter cells
  input_filtered <- input_object[,input_object[[metadata_grouping]][,1] %in% cell_types_keep]

  # Filter genes
  ligs.recs.all <- unique(c(unique(as.character(ligand_receptor_pairs$ligand)),
                     unique(as.character(ligand_receptor_pairs$receptor))))
  ligs.recs.use <- ligs.recs.all[ligs.recs.all %in% rownames(input_filtered)]
  overall.lig.rec <- Seurat::GetAssayData(input_filtered,slot="counts",assay="RNA")[ligs.recs.use,]

  ligs.keep <- apply(overall.lig.rec,1,sum)[apply(overall.lig.rec,1,sum)>min_expression &
                                              apply(overall.lig.rec,1,sum)<max_expression]

  input_filtered <- input_filtered[names(ligs.keep),]

  # Split object
  input_filtered_split <- Seurat::SplitObject(input_filtered,split=metadata_grouping)

  expr_ligs <- do.call(cbind,lapply(input_filtered_split,function(x) {

    expr <- Seurat::GetAssayData(x,assay="RNA",slot="data")
    apply(expr[rownames(expr) %in% ligand_receptor_pairs$ligand,],1,mean)

  }))

  expr_recs <- do.call(cbind,lapply(input_filtered_split,function(x) {

    expr <- Seurat::GetAssayData(x,assay="RNA",slot="data")
    apply(expr[rownames(expr) %in% ligand_receptor_pairs$receptor,],1,mean)

  }))

  return(list(expr_ligs,expr_recs))

}


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
#' @keywords internal
#' @noRd

ligand_receptor_means <- function(ligand_mean_dataframe,receptor_mean_dataframe,ligand_receptor_pairs) {

# Bind variables
Var1 <- Var2 <- interaction_pairs <- NULL

# Determine number of cell types
number_cell_types <- ncol(ligand_mean_dataframe)

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


#' Assess joint means between scrambled interaction pairs for
#' ligands and receptors across cell types
#'
#' @param ligand_receptor_means Tibble generated by ligand_receptor_means
#' function. Used to generate create a joint expression weighting for ligand
#' receptor pairs.
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
#' @param scramble_times Number of times to random pair ligands and receptors to
#' generate a background distribution.
#'
#' @return Generates a tibble containing joint resulting from scrambling ligand
#' and receptor pairs.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join
#' @importFrom stats sd
#'
#' @keywords internal
#' @noRd

ligand_receptor_scramble <- function(ligand_receptor_means_tib,ligand_mean_dataframe,receptor_mean_dataframe,ligand_receptor_pairs,scramble_times) {

# Bind variables
interaction_pairs <- scram_mean <- scram_sd <- NULL

# Determine number of cell types
number_cell_types <- ncol(ligand_mean_dataframe)

ligs.from.interact <- sapply(strsplit(colnames(ligand_receptor_means_tib)[2:ncol(ligand_receptor_means_tib)],split="_"),function(x) x[[1]])
recs.from.interact  <- sapply(strsplit(colnames(ligand_receptor_means_tib)[2:ncol(ligand_receptor_means_tib)],split="_"),function(x) x[[2]])

lig.rec.pairs.scram.mean <- lig.rec.pairs.scram.sd <- vector("list")

for (i in 1:length(ligs.from.interact)) {

	lig.use <- ligs.from.interact[i]
	rec.use <- recs.from.interact[i]

	comparisons <- expand.grid(seq(1:number_cell_types),seq(1:number_cell_types))
	scram.iter <- vector("list",length=number_cell_types)

	for (q in 1:scramble_times) {

		scramble.index <- sample(1:(number_cell_types^2),(number_cell_types^2))
		comparisons <- comparisons[scramble.index,]

		all.lig.rec <- vector("list")
		test.frame <- data.frame(ligand_mean_dataframe[lig.use,],receptor_mean_dataframe[rec.use,])

			for (z in 1:nrow(comparisons)) {

				all.lig.rec[[z]] <- mean(c(test.frame[comparisons[z,1],1],test.frame[comparisons[z,2],2]))

			}

		scram.iter[[q]] <- data.frame(do.call(rbind,all.lig.rec))

	}

	bound.scram.iter <- do.call(cbind,scram.iter)
	lig.rec.pairs.scram.mean[[i]] <- data.frame(scram_mean=apply(bound.scram.iter,1,mean))
	colnames(lig.rec.pairs.scram.mean[[i]]) <- paste(lig.use,rec.use,sep="_")
	lig.rec.pairs.scram.sd[[i]] <- data.frame(scram_sd=apply(bound.scram.iter,1,sd))
	colnames(lig.rec.pairs.scram.sd[[i]]) <- paste(lig.use,rec.use,sep="_")

}

interact.frame.scram.mean <- do.call(cbind,lig.rec.pairs.scram.mean)
interact.frame.scram.sd <- do.call(cbind,lig.rec.pairs.scram.sd)
rownames(interact.frame.scram.mean) <- rownames(interact.frame.scram.sd) <- ligand_receptor_means_tib$interaction_pairs

interact.mean.scram <- as_tibble(interact.frame.scram.mean,rownames="interaction_pairs")
interact.sd.scram <- as_tibble(interact.frame.scram.sd,rownames="interaction_pairs")

interact.mean.scram.g <- interact.mean.scram %>%
	gather("interaction","scram_mean",-interaction_pairs)
interact.sd.scram.g <- interact.sd.scram %>%
	gather("interaction","scram_sd",-interaction_pairs)

left_join(interact.mean.scram.g,interact.sd.scram.g,by=c("interaction_pairs","interaction"))

}


#' Assess statistical significance of ligand/receptor interactions within an
#' individual group of samples
#'
#' @param ligand_receptor_means Tibble generated by ligand_receptor_means
#' function. Used to generate create a joint expression weighting for ligand
#' receptor pairs.
#'
#' @param ligand_receptor_scramble Tibble generated by ligand_receptor_scramble
#' function. Used to generate scrambled values to test statistical signifiance
#' of the strength of interactions.
#'
#' @return Generates a tibble containing joint resulting from scrambling ligand
#' and receptor pairs.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#'
#' @keywords internal
#' @noRd

ligand_receptor_stats <- function(ligand_receptor_means_tib,ligand_receptor_scramble) {

# Bind variables
value <- interaction_pairs <- scram_mean <- scram_sd <- NULL

  interact.tib.g <- ligand_receptor_means_tib %>%
  	gather("interaction","value",-interaction_pairs)

  overall.lig.rec <- left_join(interact.tib.g,ligand_receptor_scramble,by=c("interaction_pairs","interaction"))

  overall.lig.rec %>%
	 mutate(p_val=1-stats::pnorm(value,mean=scram_mean,sd=scram_sd)) %>%
	 mutate(cell_type1=sapply(strsplit(interaction_pairs,split="_"),function(x) x[[1]])) %>%
	 mutate(cell_type2=sapply(strsplit(interaction_pairs,split="_"),function(x) x[[2]])) %>%
	 mutate(scram_mean=ifelse(scram_mean==0,0.1,scram_mean)) %>%
	 mutate(interact_ratio=value/scram_mean)

}
