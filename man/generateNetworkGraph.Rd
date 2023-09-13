\name{generateNetworkGraph}
\alias{generateNetworkGraph}
\alias{generateNetworkFromAdjacencyMat}
\title{
  Generate the \code{igraph} for a Network Adjacency Matrix
}
\description{
  Given the adjacency matrix of an undirected graph, returns the corresponding
  \code{\link[igraph]{igraph}} containing the list of nodes and edges.

  \code{generateNetworkFromAdjacencyMat()} is a deprecated equivalent of
  \code{generateNetworkGraph()}.
}
\usage{
generateNetworkGraph(
  adjacency_matrix
)

# Deprecated equivalent:
generateNetworkFromAdjacencyMat(
  adjacency_matrix
)
}
\arguments{
  \item{adjacency_matrix}{
    A symmetric matrix.
    Passed to
    \code{\link[igraph:graph_from_adjacency_matrix]{graph_from_adjacency_matrix()}}.
}
}
\value{
    An object of class \code{\link[igraph]{igraph}},
    containing the list of nodes and edges corresponding to \code{adjacency_matrix}.
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
toy_data <- simulateToyData(sample_size = 10)

adj_mat <-
  generateAdjacencyMatrix(
    toy_data$CloneSeq
  )

igraph <-
  generateNetworkGraph(
    adj_mat
  )
}