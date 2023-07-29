\name{getClusterStats}
\alias{getClusterStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Compute Cluster-Level Network Statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given immune repertoire sequence data and the adjacency matrix for the
  corresponding repertoire network, perform clustering and compute cluster-level
  network statistics.
}
\usage{
getClusterStats(data,
                adjacency_matrix,
                seq_col = NULL,
                count_col = NULL,
                cluster_id_col = NULL,
                degree_col = NULL,
                cluster_fun = cluster_fast_greedy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones/cells
      indexed by row.
}
  \item{adjacency_matrix}{
%%     ~~Describe \code{adjacency_matrix} here~~
    A square adjacency matrix corresponding to the repertoire network for the rows
    in \code{data}, with row dimension matching that of \code{data}. Sparse matrices
    created using the \code{Matrix} package are supported.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
      Optional. The column name or number of \code{data} containing the receptor sequences (e.g.,
      TCR CDR3 nucleotide sequence or amino acid sequence). If supplied, then related cluster-level properties will be computed.
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
      Optional. The column name or number of \code{data} containing the counts (i.e., clone count or UMI count). If supplied, then related cluster-level properties will be computed.
}
  \item{cluster_id_col}{
%%     ~~Describe \code{cluster_id_col} here~~
      The column name or number of \code{data} containing the cluster membership for each node. If \code{NULL}, clustering will be performed using \code{igraph::cluster_fast_greedy()} in order to determine cluster membership.
}
  \item{degree_col}{
%%     ~~Describe \code{degree_col} here~~
      The column name or number of \code{data} containing the network degree for each node. If \code{NULL}, the network degree will be computed using \code{igraph::degree}.
}
  \item{cluster_fun}{
%%     ~~Describe \code{net} here~~
      Passed to \code{addClusterMembership}. Controls the algorithm used for cluster identification.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The cluster-level statistics \code{diameter_length}, \code{assortativity}, \code{transitivity},
    \code{edge_density}, \code{degree_centrality_index}, \code{closeness_centrality_index},
    \code{eigen_centrality_index}, and \code{eigen_centrality_eigenvalue} are computed for each cluster
    by first generating the network graph for the cluster using the adjacency
    submatrix corresponding to that cluster. As a result, computing cluster statistics
    can potentially take some time for very large networks with many large clusters.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing one row for each cluster and the following columns:

  \item{cluster_id}{The cluster ID number.}
  \item{node_count}{The number of nodes in the cluster.}
  \item{mean_seq_length}{The mean sequence length in the cluster.}
  \item{mean_degree}{The mean network degree in the cluster.}
  \item{max_degree}{The maximum network degree in the cluster.}
  \item{seq_w_max_degree}{The receptor sequence possessing the maximum degree within the cluster.}
  \item{agg_count}{The aggregate count among all nodes in the cluster (based on the counts in \code{count_col}).}
  \item{max_count}{The maximum count among all nodes in the cluster (based on the counts in \code{count_col}).}
  \item{seq_w_max_count}{The receptor sequence possessing the maximum count within the cluster.}
  \item{diameter_length}{The longest geodesic distance in the cluster, computed using \code{length(igraph::get_diameter())}.}
  \item{assortativity}{The assortativity coefficient of the cluster's graph, based on the degree (minus one) of each node in the cluster (with the degree computed based only upon the nodes within the cluster). Computed using \code{igraph::assortativity_degree}.}
  \item{global_transitivity}{The transitivity (i.e., clustering coefficient) for the cluster's graph, which estimates the probability that adjacent vertices are connected. Computed using \code{igraph::transitivty()} with \code{type = "global"}.}
  \item{edge_density}{The number of edges in the cluster as a fraction of the maximum possible number of edges. Computed using \code{igraph::edge_density()}.}
  \item{degree_centrality_index}{The cluster-level centrality index based on degree within the cluster graph. Computed using \code{igraph::centr_degree()$centralization}.}
  \item{closeness_centrality_index}{The cluster-level centrality index based on closeness, i.e., distance to other nodes in the cluster. Computed using \code{igraph::centr_clo()$centralization}.}
  \item{eigen_centrality_index}{The cluster-level centrality index based on the eigenvector centrality scores, i.e., values of the first eigenvector of the adjacency matrix for the cluster. Computed using \code{igraph::centr_eigen()$centralization}.}
  \item{eigen_centrality_eigenvalue}{The eigenvalue corresponding to the first eigenvector of the adjacency matrix for the cluster. Computed using \code{igraph::eigen_centrality()$value}.}
}
\author{
%%  ~~who you are~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
