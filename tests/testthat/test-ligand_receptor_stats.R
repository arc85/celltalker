## Test formatting of ligand and receptor scramble
# Create dataset

filtered_lig_rec_t <- Matrix::t(filtered_lig_rec) %>% as.matrix()
lig_rec_meta <- cbind(filtered_lig_rec_t,overall_metadata)

pbmc_lig_rec <- lig_rec_meta %>%
	dplyr::filter(sample_type=="PBMC") %>%
	dplyr::select(all_of(rownames(filtered_lig_rec)),cell_types) %>%
	dplyr::group_split(cell_types)
names(pbmc_lig_rec) <- unique(lig_rec_meta$cell_types)

pbmc_lig_rec_means <- sapply(pbmc_lig_rec,function(x) {

apply(x[,1:nrow(filtered_lig_rec)],2,mean)

})

pbmc_lig_means <- pbmc_lig_rec_means[rownames(pbmc_lig_rec_means) %in% unique(ramilowski_pairs$ligand),]
pbmc_rec_means <- pbmc_lig_rec_means[rownames(pbmc_lig_rec_means) %in% unique(ramilowski_pairs$receptor),]

pbmc_joint_means <- ligand_receptor_means(ligand_mean_dataframe=pbmc_lig_means,receptor_mean_dataframe=pbmc_rec_means,ligand_receptor_pairs=ramilowski_pairs,number_cell_types=10)

pbmc_joint_scramble <- ligand_receptor_scramble(pbmc_joint_means,pbmc_lig_means,pbmc_rec_means,ramilowski_pairs,number_cell_types=10,scramble_times=10)

overall_res <- ligand_receptor_stats(pbmc_joint_means,pbmc_joint_scramble)

# Tests
test_that("structure of result", {
  expect_type(overall_res,"list")
})

test_that("dimensions", {
  expect_equal(dim(overall_res),c(6100,9))
})

test_that("expected value in row 1 column 1", {
  expect_match(overall_res[1,1] %>% pull(),"CD8 T cells_CD8 T cells")
})

test_that("value of row 6100 column 1", {
  expect_match(overall_res[6100,1] %>% pull(),"RBCs_RBCs")
})
