\name{saveNetworkPlots}
\alias{saveNetworkPlots}
\title{
    Write Plots to a PDF
}
\description{
    Given a list of plots, write all plots to a single pdf file containing one plot per page.
}
\usage{
saveNetworkPlots(
  plotlist,
  outfile = "MyRepSeqNetwork.pdf",
  pdf_width = 12,
  pdf_height = 10
)
}
\arguments{
  \item{plotlist}{
    A list whose elements can be written to pdf using the \code{\link[grDevices]{pdf}} function, such as objects of class \code{\link[ggraph]{ggraph}} or \code{\link[ggplot2]{ggplot}}.
}
  \item{outfile}{
    The file path for the pdf.
}
  \item{pdf_width}{
    Sets the page width. Passed to the \code{width} argument of the \code{\link[grDevices]{pdf}} function.
}
  \item{pdf_height}{
    Sets the page height. Passed to the \code{height} argument of the \code{\link[grDevices]{pdf}} function.
}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData()

net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )

net$plots <-
  generateNetworkGraphPlots(
    net$igraph,
    net$node_data,
    color_nodes_by =
      c("SampleID", "CloneCount")
  )

saveNetworkPlots(
  net$plots,
  outfile =
    file.path(tempdir(), "network.pdf")
)


# clean up temp directory
unlink(list.dirs(tempdir(), full.names = TRUE),
       recursive = TRUE)
}