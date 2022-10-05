#' Creates a circos plot from the list of ligands and receptors
#'
#' @param ligand_receptor_frame Resulting tibble (usually filtered in some way)
#' from the celltalk function.
#'
#' @param cell_group_colors Colors used for the groups of cells in the outer
#' track of the circos plot.
#'
#' @param ligand_color Color to use for ligands. Defaults to "blue".
#'
#' @param receptor_color Color to use for the receptors. Defaults to "red".
#'
#' @param cex_outer Size of the text for the cell groups in the outer layer of
#' the circos plot. Default is 0.5.
#'
#' @param cex_inner Size of the text for the ligand and receptors in the
#' inner layer of the circos plot. Default is 0.4.
#'
#' @return Generates a circos plot connecting ligands and receptors across cell types for a given sample group
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_split
#' @importFrom dplyr summarize
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom circlize CELL_META
#' @importFrom circlize circos.link
#' @importFrom circlize circos.track
#' @importFrom circlize circos.clear
#' @importFrom circlize circos.par
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.rect
#' @importFrom circlize circos.text
#'
#' @export

circos_plot <- function(ligand_receptor_frame,
  cell_group_colors,
  ligand_color="blue",
  receptor_color="red",
  cex_outer=0.5,
  cex_inner=0.4,
  inner_facing="downward",
  inner_niceFacing=FALSE) {

  # Bind variables
  cell_type1 <- lig <- cell_type2 <- rec <- classes <- ranges <-
    max_range <- to_class <- to_rec <- lig_rec <- ordered_lig_rec <- type <-
    lig.rec <- to.class <- to.rec <- ordered.lig.rec <- NULL

  # Reformat data
  part1 <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type1,lig) %>%
    mutate(type="lig")
  part2 <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type2,rec) %>%
    mutate(type="rec")
  colnames(part1) <- colnames(part2) <- c("classes","lig.rec","type")

  part12 <- rbind(part1,part2) %>%
    group_by(classes) %>%
    group_split()

  part12 <- lapply(part12,function(x) {

    x <- x %>%
      mutate(ordered.lig.rec=paste(type,lig.rec,sep="_")) %>%
      mutate(ranges=as.numeric(as.factor(ordered.lig.rec))) %>%
      select(-ordered.lig.rec)

  })

  part12 <- do.call(rbind,part12)

  to.join <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type1,cell_type2,lig,rec)

  colnames(to.join)[1:3] <- c("classes","to.class","lig.rec")

  part3 <- part12

  joined <- left_join(part3,to.join,by=c("classes","lig.rec"))
  joined$to.rec <- NA

  for (i in 1:nrow(joined)) {

    sub.group <- joined[i,]
    sub.joined <- joined %>% filter(classes==sub.group$to.class)
    joined$to.rec[i] <- sub.joined[match(sub.group$rec,sub.joined$lig.rec),"ranges"] %>% pull()

  }

  final.construct <- joined

  # Repair single class
  single.class <- final.construct %>%
    group_by(classes) %>%
    summarize(max_range=max(ranges)) %>%
    filter(max_range==1) %>%
    pull(classes)

  if (!length(single.class)==0) {

    for (i in 1:length(single.class)) {

      row.add <- final.construct[final.construct$classes==single.class[i],][1,]
      row.add$ranges <- 2

      final.construct <- rbind(final.construct,row.add)
      final.construct <- final.construct %>%
        arrange(classes)

    }

  }

  final.construct <- final.construct %>%
    arrange(classes,ranges)

  circos.clear()
  circos.par(gap.degree=10,track.margin=c(0,0.2))
  circos.initialize(factors=final.construct$classes,x=final.construct$ranges)

  suppressMessages({
  circos.track(ylim = c(0, 1),track.height=0.1,panel.fun = function(x, y) {
    circos.rect(CELL_META$cell.xlim[1],CELL_META$cell.ylim[1],CELL_META$cell.xlim[2],CELL_META$cell.ylim[2],col=cell_group_colors[CELL_META$sector.numeric.index])
    circos.text(CELL_META$xcenter, y=2.5, CELL_META$sector.index,
                facing = "downward",cex=cex_outer)
  })
  })

  ## Build interior track with ligand/receptors colors and gene labels
  circos.track(ylim = c(0, 1),track.height=0.05,bg.border="white")

  # Define multiplers for each sector
  final.construct2 <- final.construct %>%
    select(classes,lig.rec,ranges,type) %>%
    distinct() %>%
    arrange(classes,ranges)
  ref.tab <- unname(table(final.construct2$classes))
  sec.multi <- (ref.tab-1)/ref.tab
  names(sec.multi) <- names(table(final.construct2$classes))

  # Loop to construct all sectors
  # Ligands first
  # Split into list of sectors
  int.types.list <- final.construct2 %>%
    group_split(classes)

  names(int.types.list) <- sapply(int.types.list,function(x) x$classes[1])

  int.types.list.multi <- int.types.list.individ <- list("NA")

  int.types.list.multi <- int.types.list[!names(int.types.list) %in% single.class]

  int.types.list.individ <- int.types.list[names(int.types.list) %in% single.class]

  for (i in 1:length(int.types.list.multi)) {

    for (a in 1:nrow(int.types.list.multi[[i]])) {

      if (a==1) {

        sec.multi.use <- sec.multi[names(sec.multi)==int.types.list.multi[[i]]$classes[1]]

        suppressMessages({
        circos.rect(1,0,1+sec.multi.use*a,1,sector.index=int.types.list.multi[[i]]$classes[a],
                    col=ifelse(int.types.list.multi[[i]]$type[a]=="lig",ligand_color,receptor_color),track.index = 2)
        circos.text(1+sec.multi.use*a/2,4,sector.index=int.types.list.multi[[i]]$classes[a],
                    labels=int.types.list.multi[[i]]$lig.rec[a],track.index = 2,facing=inner_facing,niceFacing=inner_niceFacing,cex=cex_inner)
        })
      } else {

        sec.multi.use <- sec.multi[names(sec.multi)==int.types.list.multi[[i]]$classes[1]]

        suppressMessages({
        circos.rect(1+sec.multi.use*(a-1),0,1+sec.multi.use*a,1,sector.index=int.types.list.multi[[i]]$classes[a],
                    col=ifelse(int.types.list.multi[[i]]$type[a]=="lig",ligand_color,receptor_color),track.index = 2)
        circos.text(1+sec.multi.use*a-sec.multi.use/2,4,sector.index=int.types.list.multi[[i]]$classes[a],
                    labels=int.types.list.multi[[i]]$lig.rec[a],track.index = 2,facing=inner_facing,niceFacing=inner_niceFacing,cex=cex_inner)
        })
      }

    }

  }

  if (length(int.types.list.individ)>0) {

    for (i in 1:length(int.types.list.individ)) {

      circos.rect(1,0,2,1,sector.index=int.types.list.individ[[i]]$classes[1],
                  col=ifelse(int.types.list.individ[[i]]$type[1]=="lig",ligand_color,receptor_color),track.index = 2)
      circos.text(1.5,4,sector.index=int.types.list.individ[[i]]$classes[1],
                  labels=int.types.list.individ[[i]]$lig.rec[1],track.index = 2,facing=inner_facing,niceFacing=inner_niceFacing,cex=cex_inner)

    }

  }


  ## Draw links
  final.construct3 <- final.construct %>%
    select(classes,lig.rec,ranges,to.class,to.rec) %>%
    distinct()

  int.types.list <- final.construct3 %>%
    group_split(classes)

  names(int.types.list) <- sapply(int.types.list,function(x) x$classes[1])


  for (i in 1:length(int.types.list)) {

    for (a in 1:nrow(int.types.list[[i]])) {

      target <- which(!is.na(match(names(int.types.list),int.types.list[[i]]$to.class[[a]])))

      if (length(target)==0) {

      } else {

        if (!int.types.list[[i]]$to.class[[a]] %in% single.class) {

          circos.link(int.types.list[[i]]$classes[a], 1+sec.multi[i]*int.types.list[[i]]$ranges[a]-sec.multi[i]/2,
                      int.types.list[[i]]$to.class[[a]], 1+sec.multi[target]*int.types.list[[i]]$to.rec[a]-sec.multi[target]/2,
                      0.43, 0.43, directional=1, lwd=3, arr.length=0.2, arr.width=(3*0.1)/2)

        } else {

          circos.link(int.types.list[[i]]$classes[a], 1+sec.multi[i]*int.types.list[[i]]$ranges[a]-sec.multi[i]/2,
                      int.types.list[[i]]$to.class[[a]], 1.5,
                      0.43, 0.43, directional=1, lwd=3, arr.length=0.2, arr.width=(3*0.1)/2)
        }


      }

    }

  }

}
