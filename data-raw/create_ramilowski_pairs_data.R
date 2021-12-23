## Code used to prepare the `ramilowski_pairs.rda` dataset goes here

# Read in curated Supplementary Table 3 from "A draft network of ligand-receptor-mediated multicellular signalling in human" from Nature Communications, July 2015 by Jordan Ramilowski et al. This table contains ligands, receptors, and receptor/ligand pairs

ramilowski_pairs <- read.csv("ramilowski-supplementary_s3_all_pairs.txt",sep="\t")
ramilowski_pairs <- ramilowski_pairs[,c(2,3,1)]
colnames(ramilowski_pairs) <- c("ligand","receptor","pair")

usethis::use_data(ramilowski_pairs)
