\name{addNodeNetworkStats}
\alias{addNodeNetworkStats}
\title{
  Compute Node-Level Network Properties
}
\description{
  Compute a specified set of node-level network properties for a given network. Each property is added as a variable to the node-level metadata for the network.

}
\usage{
addNodeNetworkStats(
  data,
  net,
  stats_to_include = chooseNodeStats(),
  cluster_fun = cluster_fast_greedy
)
}
\arguments{
  \item{data}{
      A data frame containing the node-level metadata for the network, such as the data frame named \code{node_data} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{net}{
      The \code{\link[igraph]{igraph}} object for the network, such as the element named \code{igraph} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{stats_to_include}{
      Specifies which network properties to compute. Either a vector created using \code{\link{chooseNodeStats}} or \code{\link{exclusiveNodeStats}}, or \code{"all"} to compute all network properties, or \code{"cluster_id_only"} to compute only the \code{cluster_id} property.
}
  \item{cluster_fun}{
      The function to use for the clustering algorithm. Passed to \code{\link{addClusterMembership}}. See \code{\link{clustering_algorithms}}.
}
}
\details{
  Node-level network properties are properties that pertain to each individual node in the network graph.

  Some are local properties, meaning that their value for a given node depends only on a subset of the nodes in the network. One example is the network degree of a given node, which represents the number of other nodes that are directly joined to the given node by an edge connection.

  Other properties are global properties, meaning that their value for a given node depends on all of the nodes in the network. An example is the authority score of a node, which is computed using the entire graph adjacency matrix (if we denote this matrix by \eqn{A}, then the principal eigenvector of \eqn{A^T A} represents the authority scores of the network nodes).

  See \code{\link{chooseNodeStats}} for a list of available node-level network properties.
}
\value{
  A copy of \code{data} with with an additional column for each network property computed. See \code{\link{chooseNodeStats}} for the network property names, which are used as the column names.
}
\seealso{
\code{\link{chooseNodeStats}}
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
toy_data <- simulateToyData()

# Generate network
net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )

# Add default set of node properties
data_w_default_stats <-
  addNodeNetworkStats(
    net$node_data, net$igraph
  )

# Modify default set of node properties
data_w_custom_stats1 <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include =
      chooseNodeStats(
        closeness = TRUE,
        page_rank = FALSE
      )
  )

# Add only the spepcified node properties
data_w_custom_stats2 <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include =
      exclusiveNodeStats(
        degree = TRUE,
        transitivity = TRUE
      )
  )

# Add all node-level network properties
data_w_all_stats <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include = "all"
  )

}
}