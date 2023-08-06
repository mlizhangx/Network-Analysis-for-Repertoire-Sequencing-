\name{saveNetwork}
\alias{saveNetwork}
\title{
    Save List of Network Objects
}
\description{
    Given a list of network objects such as that returned by \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}, saves its contents according to the specified file format scheme.
}
\usage{
saveNetwork(
  net,
  output_dir = getwd(),
  output_type = "individual",
  output_filename = "MyRepSeqNetwork",
  pdf_width = 12,
  pdf_height = 10
)
}
\arguments{
  \item{net}{
    A list conforming to the format returned by \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}. See details.
}
  \item{output_dir}{
    A file path specifying the directory in which to write the file(s).
}
  \item{output_type}{
  A character string specifying the file format scheme to use when writing output to file. Valid options are \code{"individual"}, \code{"rds"} and \code{"rda"}. See detials.
}
  \item{output_filename}{
  A character string. All files saved will have file names beginning with this value.
}
  \item{pdf_width}{
    If the list contains plots, this controls the width of each plot when writing to pdf. Passed to the \code{width} argument of the \code{\link[grDevices]{pdf}} function.
}
  \item{pdf_height}{
    If the list contains plots, this controls the height of each plot when writing to pdf. Passed to the \code{height} argument of the \code{\link[grDevices]{pdf}} function.
}
}
\details{
    The list provided to the \code{net} argument must contain at least the elements \code{igraph} (of class \code{\link[igraph]{igraph}}), \code{adjacency_matrix} (of class \code{matrix} or \code{dgCMatrix}), and \code{node_data} (a data frame). The components \code{cluster_data} (a data frame) and \code{plots} (a list of \code{\link[ggraph]{ggraph}} objects) are optional, and will also be saved if present.

    By default, each list element is saved as a separate, uncompressed file: data frames, including \code{node_data} and \code{cluster_data}, are saved as .csv files; \code{igraph} is saved as a .txt (using \code{\link[igraph]{write_graph}} with \code{format = "edgelist"}); \code{adjacency_matrix} is saved as a .mtx file (for the \code{dgCMatrix} class) or .csv file (for the \code{Matrix} class); the elements of \code{plots} are printed to a .pdf file using \code{\link[grDevices]{pdf}}, with the pdf file containing one plot per page.

    Setting \code{output_type} to \code{"rda"} or \code{"rds"} writes the entire  list \code{net} to a single compressed file (a .RData file or a .RDS file, respectively), while also printing the plots to a .pdf file.

}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
\dontrun{
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")

# Save network objects
saveNetwork(net, output_dir = tempdir())
}
}