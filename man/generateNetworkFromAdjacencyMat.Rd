\name{generateNetworkFromAdjacencyMat}
\alias{generateNetworkFromAdjacencyMat}
\title{
  Generate an igraph From a Network Adjacency Matrix
}
\description{
  Given the adjacency matrix of an undirected graph, returns the corresponding \code{\link[igraph]{igraph}} object containing the list of nodes and edges.
}
\usage{
generateNetworkFromAdjacencyMat(adjacency_matrix)
}
\arguments{
  \item{adjacency_matrix}{
    A symmetric matrix. Passed to \code{\link[igraph]{graph_from_adjacency_matrix}}.
}
}
\value{
    An object of class \code{\link[igraph]{igraph}}.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData(sample_size = 10)

adjmat <-
  sparseAdjacencyMatFromSeqs(
    toy_data$CloneSeq
  )

net <-
  generateNetworkFromAdjacencyMat(
    adjmat
  )
}