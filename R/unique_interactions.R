#' Unique interactions between two groups
#'
#' Identifies interaction pairs that are exclusively present in one group versus another at the group level
#'
#' @param putative.interactions.tib A tibble output from \code{putative_interactions}.
#' @param group1 Sample group to evaluate for unqiue ligands or receptors participating in interactions.
#' @param group2 Sample group as a comparitor to group1
#' @param ligands.and.receptors A matrix with 3 columns, the first of which corresponds to ligands, the second of which corresponds to receptors and the last of which corresponds to ligand_receptor interactions. Column names should be "ligand", "receptor", and "pairs" respectively.
#'
#' @return Returns a tibble containing interactions unique to group1 versus group2, group2 versus group1 and interactions that are common between the groups.
#'
#' @export

unique_interactions <- function(putative.interactions.tib,group1,group2,ligands.and.receptors) {

group1.exp <- unnest(putative.interactions.tib[putative.interactions.tib$group %in% group1,2]) %>% pull()
group2.exp <- unnest(putative.interactions.tib[putative.interactions.tib$group %in% group2,2]) %>% pull()
names(group1.exp) <- names(group2.exp) <- ligands.and.receptors$pair

non.null1 <- group1.exp[sapply(group1.exp,function(x) !any(is.null(x[[1]])))]
non.null2 <- group2.exp[sapply(group2.exp,function(x) !any(is.null(x[[1]])))]

unique1 <- setdiff(names(non.null1),names(non.null2))
unique2 <- setdiff(names(non.null2),names(non.null1))

common.interactions <- intersect(names(non.null1),names(non.null2))

enframe(list("unique1v2"=unique1,"unique2v1"=unique2,"common"=common.interactions),name="comparison",value="ligands.and.receptors")

}
