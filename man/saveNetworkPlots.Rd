\name{saveNetworkPlots}
\alias{saveNetworkPlots}
\title{
    Write Plots to a PDF
}
\description{
    Given a list of plots, write all plots to a single pdf file containing one
    plot per page, and optionally save the graph layout as a csv file.
}
\usage{
saveNetworkPlots(
  plotlist,
  outfile,
  pdf_width = 12,
  pdf_height = 10,
  outfile_layout = NULL,
  verbose = FALSE
)
}
\arguments{
  \item{plotlist}{
    A named list whose elements are of class \code{\link[ggraph]{ggraph}}. May
    also contain an element named \code{graph_layout} with the matrix specifying
    the graph layout.
}
  \item{outfile}{
    A \code{\link{connection}} or a character string containing the file path
    used to save the pdf.
}
  \item{pdf_width}{
    Sets the page width. Passed to the \code{width} argument of
    \code{\link[grDevices:pdf]{pdf()}}.
}
  \item{pdf_height}{
    Sets the page height. Passed to the \code{height} argument of
    \code{\link[grDevices:pdf]{pdf()}}.
}
  \item{outfile_layout}{
    An optional \code{\link{connection}} or file path for saving the graph layout.
    Passed to
    the \code{file} argument of \code{\link[=write]{write()}},
    which is called with \code{ncolumns = 2}.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
}
\value{
Returns \code{TRUE}, invisibly.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
set.seed(42)
toy_data <- simulateToyData()

net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )

net <-
  addPlots(
    net,
    color_nodes_by =
      c("SampleID", "CloneCount"),
    print_plots = TRUE
  )

saveNetworkPlots(
  net$plots,
  outfile =
    file.path(tempdir(), "network.pdf"),
  outfile_layout =
    file.path(tempdir(), "graph_layout.txt")
)

# Load saved graph layout
graph_layout <- matrix(
  scan(file.path(tempdir(), "graph_layout.txt"), quiet = TRUE),
  ncol = 2
)
all.equal(graph_layout, net$plots$graph_layout)


\dontshow{
# clean up temporary directory
file.remove(
  file.path(tempdir(), c("network.pdf", "graph_layout.txt"))
)
}
}