\name{getClusterStats}
\alias{getClusterStats}
\title{
    Compute Cluster-Level Network Statistics
}
\description{
  Given the node-level metadata and adjacency matrix for an immune repertoire network, perform clustering and compute cluster-level network properties.
}
\usage{
getClusterStats(
    data,
    adjacency_matrix,
    seq_col = NULL,
    count_col = NULL,
    cluster_id_col = NULL,
    degree_col = NULL,
    cluster_fun = cluster_fast_greedy
)
}
\arguments{
  \item{data}{
      A data frame containing the node-level metadata for the network, with each row corresponding to a network node.
}
  \item{adjacency_matrix}{
    The adjacency matrix for the network.
}
  \item{seq_col}{
    Specifies the column of \code{data} containing the receptor sequences. Accepts either the column name or column index. If provided, then related cluster-level properties will be computed.
}
  \item{count_col}{
      Specifies the column of \code{data} containing a variable with counts (such as clone count or UMI count). Accepts either the column name or column index. If provided, then related cluster-level properties will be computed.
}
  \item{cluster_id_col}{
      If cluster analysis has been previously performed to partition the graph into clusters, this argument specifies the column of \code{data} containing the ID of the cluster to which each node belongs. Accepts either the column name or column index. If \code{NULL}, then \code{\link{addClusterMembership}} will be used to partition the graph into clusters.
}
  \item{degree_col}{
      Specifies the column of \code{data} containing the network degree of each node. Accepts either the column name or column index. If \code{NULL}, the network degree of each node will be computed.
}
  \item{cluster_fun}{
      If \code{cluster_id_col = NULL}, this argument specifies the clustering algorithm used to partition the graph into clusters. Passed to \code{\link{addClusterMembership}}.
}
}
\details{
    The cluster-level properties \code{diameter_length}, \code{assortativity}, \code{transitivity}, \code{edge_density}, \code{degree_centrality_index}, \code{closeness_centrality_index}, \code{eigen_centrality_index}, and \code{eigen_centrality_eigenvalue} are computed for each cluster using the corresponding submatrix of the network adjacency matrix.


}
\value{
    A data frame containing one row for each cluster and the following variables:

  \item{cluster_id}{The cluster ID number.}
  \item{node_count}{The number of nodes in the cluster.}
  \item{mean_seq_length}{The mean sequence length in the cluster.}
  \item{mean_degree}{The mean network degree in the cluster.}
  \item{max_degree}{The maximum network degree in the cluster.}
  \item{seq_w_max_degree}{The receptor sequence possessing the maximum degree within the cluster.}
  \item{agg_count}{The aggregate count among all nodes in the cluster (based on the counts in \code{count_col}).}
  \item{max_count}{The maximum count among all nodes in the cluster (based on the counts in \code{count_col}).}
  \item{seq_w_max_count}{The receptor sequence possessing the maximum count within the cluster.}
  \item{diameter_length}{The longest geodesic distance in the cluster, computed as the length of the vector returned by \code{\link[igraph]{get_diameter}}.}
  \item{assortativity}{The assortativity coefficient of the cluster's graph, based on the degree (minus one) of each node in the cluster (with the degree computed based only upon the nodes within the cluster). Computed using \code{\link[igraph]{assortativity_degree}}.}
  \item{global_transitivity}{The transitivity (i.e., clustering coefficient) for the cluster's graph, which estimates the probability that adjacent vertices are connected. Computed using \code{\link[igraph]{transitivity}} with \code{type = "global"}.}
  \item{edge_density}{The number of edges in the cluster as a fraction of the maximum possible number of edges. Computed using \code{\link[igraph]{edge_density}}.}
  \item{degree_centrality_index}{The cluster-level centrality index based on degree within the cluster graph. Computed as the \code{centralization} element of the output from  \code{\link[igraph]{centr_degree}}.}
  \item{closeness_centrality_index}{The cluster-level centrality index based on closeness, i.e., distance to other nodes in the cluster. Computed using \code{\link[igraph]{centralization}}.}
  \item{eigen_centrality_index}{The cluster-level centrality index based on the eigenvector centrality scores, i.e., values of the first eigenvector of the adjacency matrix for the cluster. Computed as the \code{centralization} element of the output from  \code{\link[igraph]{centr_eigen}}.}
  \item{eigen_centrality_eigenvalue}{The eigenvalue corresponding to the first eigenvector of the adjacency matrix for the cluster. Computed as the \code{value} element of the output from  \code{\link[igraph]{eigen_centrality}}.}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

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

# Get cluster stats
cluster_info <-
  getClusterStats(net$node_data, net$adjacency_matrix,
                  seq_col = "CloneSeq", count_col = "CloneCount")
}
}
