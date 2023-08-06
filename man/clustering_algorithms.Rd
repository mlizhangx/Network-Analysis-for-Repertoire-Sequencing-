\name{Clustering Algorithms}
\alias{clustering_algorithms}
\alias{cluster_edge_betweenness}
\alias{cluster_fast_greedy}
\alias{cluster_fluid_communities}
\alias{cluster_infomap}
\alias{cluster_label_prop}
\alias{cluster_leading_eigen}
\alias{cluster_leiden}
\alias{cluster_louvain}
\alias{cluster_optimal}
\alias{cluster_spinglass}
\alias{cluster_walktrap}
\title{Available Clustering Algorithms}
\description{
Clustering algorithm functions available for use within the \code{NAIR} package.
}
\details{
The available functions that can be used as arguments to the \code{cluster_fun} argument of package functions such as \code{\link{buildRepSeqNetwork}} and \code{\link{addClusterMembership}} are listed below:

\code{\link[igraph]{cluster_edge_betweenness}}

\code{\link[igraph]{cluster_fast_greedy}}

\code{\link[igraph]{cluster_fluid_communities}}

\code{\link[igraph]{cluster_infomap}}

\code{\link[igraph]{cluster_label_prop}}

\code{\link[igraph]{cluster_leading_eigen}}

\code{\link[igraph]{cluster_leiden}}

\code{\link[igraph]{cluster_louvain}}

\code{\link[igraph]{cluster_optimal}}

\code{\link[igraph]{cluster_spinglass}}

\code{\link[igraph]{cluster_walktrap}}

These functions have been imported from the \code{\link[igraph]{igraph}} package. Click on one of the functions to see its original documentation.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}