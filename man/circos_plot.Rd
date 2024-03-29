% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circos_plot.R
\name{circos_plot}
\alias{circos_plot}
\title{Creates a circos plot from the list of ligands and receptors}
\usage{
circos_plot(
  ligand_receptor_frame,
  cell_group_colors,
  ligand_color = "blue",
  receptor_color = "red",
  cex_outer = 0.5,
  cex_inner = 0.4
)
}
\arguments{
\item{ligand_receptor_frame}{Resulting tibble (usually filtered in some way)
from the celltalk function.}

\item{cell_group_colors}{Colors used for the groups of cells in the outer
track of the circos plot.}

\item{ligand_color}{Color to use for ligands. Defaults to "blue".}

\item{receptor_color}{Color to use for the receptors. Defaults to "red".}

\item{cex_outer}{Size of the text for the cell groups in the outer layer of
the circos plot. Default is 0.5.}

\item{cex_inner}{Size of the text for the ligand and receptors in the
inner layer of the circos plot. Default is 0.4.}
}
\value{
Generates a circos plot connecting ligands and receptors across cell types for a given sample group
}
\description{
Creates a circos plot from the list of ligands and receptors
}
