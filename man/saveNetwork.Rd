\name{saveNetwork}
\alias{saveNetwork}
\title{
    Save List of Network Objects
}
\description{
    Given a list of network objects such as that returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or \code{\link{generateNetworkObjects}},
    saves its contents according to the specified file format scheme.
}
\usage{
saveNetwork(
  net,
  output_dir,
  output_type = "rds",
  output_name = "MyRepSeqNetwork",
  pdf_width = 12,
  pdf_height = 10,
  verbose = FALSE,
  output_filename = deprecated()
)
}
\arguments{
  \item{net}{
      A list of network objects returned by
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
}
  \item{output_dir}{
    A file path specifying the directory in which to write the file(s).
}
  \item{output_type}{
  A character string specifying the file format scheme to use when writing output
  to file. Valid options are \code{"individual"}, \code{"rds"} and \code{"rda"}.
  See detials.
}
  \item{output_name}{
  A character string. All files saved will have file names beginning with this value.
}
  \item{pdf_width}{
    If the list contains plots, this controls the width of each plot when writing
    to pdf. Passed to the \code{width} argument of the \code{\link[grDevices]{pdf}}
    function.
}
  \item{pdf_height}{
    If the list contains plots, this controls the height of each plot when writing
    to pdf. Passed to the \code{height} argument of the \code{\link[grDevices]{pdf}}
    function.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
  \item{output_filename}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
  Equivalent to \code{output_name}.
}
}
\details{
    The list \code{net} must contain the named elements
    \code{igraph} (of class \code{\link[igraph]{igraph}}),
    \code{adjacency_matrix} (a \code{\link{matrix}} or
    \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} encoding edge connections),
    and \code{node_data} (a \code{\link{data.frame}} containing node metadata),
    all corresponding to the same network. The list returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} and
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}
    is an example of a valid input for the \code{net} argument.

    The additional elements \code{cluster_data} (a \code{\link{data.frame}}) and
    \code{plots} (a list containing objects of class
    \code{\link[ggraph:ggraph]{ggraph}}
    and possibly one \code{\link{matrix}} named \code{graph_layout})
    will also be saved, if present.

    By default, the list \code{net} is saved to a compressed data file in the
    RDS format, while any plots present are printed to a single pdf containing
    one plot per page.

    The name of each saved file begins with the value of \code{output_name}.
    When \code{output_type} is one of \code{"rds"} or \code{"rda"},
    only two files are saved (the rds/rda and the pdf); for each file,
    \code{output_name} is followed by the appropriate file extension.

    When \code{output_type = "individual"}, each element of \code{net} is saved
    as a separate file, where \code{output_name} is followed by:
    \itemize{
    \item \code{_NodeMetadata.csv} for \code{node_data}
    \item \code{_ClusterMetadata.csv} for \code{cluster_data}
    \item \code{_EdgeList.txt} for \code{igraph}
    \item \code{_AdjacencyMatrix.mtx} for \code{adjacency_matrix}
    \item \code{_Plots.rda} for \code{plots}
    \item \code{_GraphLayout.txt} for \code{plots$graph_layout}
    \item \code{_Details.rds} for \code{details}
    }

    \code{node_data} and \code{cluster_data} are saved using
    \code{\link[utils:write.csv]{write.csv()}},
    with \code{row.names} being \code{TRUE} for \code{node_data}
    and \code{FALSE} for \code{cluster_data}.
    The \code{igraph} is saved using
    \code{\link[igraph:write_graph]{write_graph()}} with \code{format = "edgelist"}.
    The adjacency matrix is saved using \code{\link[Matrix:writeMM]{writeMM()}}.
    The graph layout is saved using \code{\link[=write]{write()}} with
    \code{ncolumns = 2}.
}
\value{
Returns \code{TRUE} if output is saved, otherwise returns \code{FALSE} (with a warning if \code{output_dir} is non-null and the specified directory does not exist and could not be created).
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

net <- buildRepSeqNetwork(
  toy_data,
  seq_col = "CloneSeq",
  node_stats = TRUE,
  cluster_stats = TRUE,
  color_nodes_by = c("transitivity", "SampleID")
)

# save as single RDS file
saveNetwork(
  net,
  output_dir = tempdir(),
  verbose = TRUE
)

saveNetwork(
  net,
  output_dir = tempdir(),
  output_type = "individual",
  verbose = TRUE
)

\dontshow{
# clean up temp files
file.remove(
  file.path(
    tempdir(),
    c("MyRepSeqNetwork_NodeMetadata.csv",
      "MyRepSeqNetwork_ClusterMetadata.csv",
      "MyRepSeqNetwork_EdgeList.txt",
      "MyRepSeqNetwork_AdjacencyMatrix.mtx",
      "MyRepSeqNetwork_Details.rds",
      "MyRepSeqNetwork_Plots.rda",
      "MyRepSeqNetwork_GraphLayout.txt",
      "MyRepSeqNetwork.pdf",
      "MyRepSeqNetwork.rds"
    )
  )
)
}
}