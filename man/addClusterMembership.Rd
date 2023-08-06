\name{addClusterMembership}
\alias{addClusterMembership}
\title{
    Perform Cluster Analysis on the Nodes of a Network Graph
}
\description{
    Given the \code{igraph} and node-level metadata for a network, partitions the graph into densely-connected subgraphs (clusters) using the specified clustering algorithm, and adds a cluster membership ID variable to the node-level metadata.
}
\usage{
addClusterMembership(
  data,
  net,
  fun = cluster_fast_greedy
)
}
\arguments{
  \item{data}{
      A data frame containing the node-level metadata for the network, such as the data frame named \code{node_data} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{net}{
      The \code{\link[igraph]{igraph}} object for the network, such as the element named \code{igraph} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{fun}{
      The function to use for the clustering algorithm. See \code{\link{clustering_algorithms}}.
}
}
\details{
    Cluster analysis involves using a community-finding algorithm to partition the network graph into clusters (densely-connected subgraphs). These clusters represent groups of clones/cells with similar receptor sequences.

    Each cluster is assigned a numeric cluster ID, and the cluster membership ID of each node is recorded. The cluster membership ID is the cluster ID of the cluster to which the node belongs.

    Refer to the help topic \code{\link{clustering_algorithms}} for information on the available clustering algorithms that can be used to perform cluster analysis.
}
\value{
    A copy of \code{data} with an additional column named \code{cluster_id} that contains the cluster membership ID for the network node corresponding to each row.
}
\seealso{
\code{\link{addClusterLabels}}
\code{\link{getClusterStats}}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData()

# Generate network
net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )

# Perform cluster analysis,
# add cluster membership ID
data_w_clusterID <-
  addClusterMembership(
    net$node_data,
    net$igraph
  )
}
