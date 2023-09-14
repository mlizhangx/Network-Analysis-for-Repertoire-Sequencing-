\name{getClusterStats}
\alias{getClusterStats}
\title{
    Compute Cluster-Level Network Properties
}
\description{
  Given the node-level metadata and adjacency matrix for a network graph that has
  been partitioned into clusters, computes network properties for the clusters
  and returns them in a data frame.

  \code{\link[=addClusterStats]{addClusterStats()}}
  is preferred to \code{getClusterStats()} in most situations.
}
\usage{
getClusterStats(
  data,
  adjacency_matrix,
  seq_col = NULL,
  count_col = NULL,
  cluster_id_col = "cluster_id",
  degree_col = NULL,
  cluster_fun = deprecated(),
  verbose = FALSE
)
}
\arguments{
  \item{data}{
      A data frame containing the node-level metadata for the network, with each
      row corresponding to a network node.
}
  \item{adjacency_matrix}{
    The adjacency matrix for the network.
}
  \item{seq_col}{
      Specifies the column(s) of \code{data} containing
      the receptor sequences upon whose similarity the network is based.
      Accepts a character or numeric vector of length 1 or 2,
      containing either column names or column indices.
      If provided, then related cluster-level properties will be computed.
}
  \item{count_col}{
      Specifies the column of \code{data} containing
      a measure of abundance (such as clone count or UMI count).
      Accepts a character string containing the column name
      or a numeric scalar containing the column index.
      If provided, related cluster-level properties will be computed.
}
  \item{cluster_id_col}{
      Specifies the column of \code{data} containing the cluster membership
      variable that identifies the cluster to which each node belongs.
      Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
  \item{degree_col}{
      Specifies the column of \code{data} containing
      the network degree of each node.
      Accepts a character string containing the column name
      or a numeric scalar containing the column index.
      If the column does not exist, the network degree will be computed.
}
  \item{cluster_fun}{
      \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
      Does nothing.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
}
\details{
    To use \code{getClusterStats()},
    the network graph must first be partitioned into clusters,
    which can be done using
    \code{\link[=addClusterMembership]{addClusterMembership()}}.
    The name of the cluster membership variable in the node metadata
    must be provided to the \code{cluster_id_col} argument
    when calling \code{getClusterStats()}.
}
\value{
    A data frame containing one row for each cluster in the network
    and the following variables:

  \item{cluster_id}{
  The cluster ID number.
  }
  \item{node_count}{
  The number of nodes in the cluster.
  }
  \item{mean_seq_length}{
  The mean sequence length in the cluster.
  Only present when \code{length(seq_col) == 1}.
  }
  \item{A_mean_seq_length}{
  The mean first sequence length in the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{B_mean_seq_length}{
  The mean second sequence length in the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{mean_degree}{
  The mean network degree in the cluster.
  }
  \item{max_degree}{
  The maximum network degree in the cluster.
  }
  \item{seq_w_max_degree}{
  The receptor sequence possessing the maximum degree within the cluster.
  Only present when \code{length(seq_col) == 1}.
  }
  \item{A_seq_w_max_degree}{
  The first sequence of the node possessing the maximum degree within the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{B_seq_w_max_degree}{
  The second sequence of the node possessing the maximum degree within the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{agg_count}{
  The aggregate count among all nodes in the cluster (based on the counts in
  \code{count_col}).
  }
  \item{max_count}{
  The maximum count among all nodes in the cluster (based on the counts in
  \code{count_col}).
  }
  \item{seq_w_max_count}{
  The receptor sequence possessing the maximum count within the cluster.
  Only present when \code{length(seq_col) == 1}.
  }
  \item{A_seq_w_max_count}{
  The first sequence of the node possessing the maximum count within the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{B_seq_w_max_count}{
  The second sequence of the node possessing the maximum count within the cluster.
  Only present when \code{length(seq_col) == 2}.
  }
  \item{diameter_length}{
  The longest geodesic distance in the cluster, computed as the length of the
  vector returned by \code{\link[igraph:get_diameter]{get_diameter()}}.
  }
  \item{assortativity}{
  The assortativity coefficient of the cluster's graph, based on the degree
  (minus one) of each node in the cluster (with the degree computed based only
  upon the nodes within the cluster). Computed using
  \code{\link[igraph:assortativity_degree]{assortativity_degree()}}.
  }
  \item{global_transitivity}{
  The transitivity (i.e., clustering coefficient) for the cluster's graph, which
  estimates the probability that adjacent vertices are connected. Computed using
  \code{\link[igraph:transitivity]{transitivity()}} with \code{type = "global"}.
  }
  \item{edge_density}{
  The number of edges in the cluster as a fraction of the maximum possible number
  of edges. Computed using \code{\link[igraph:edge_density]{edge_density()}}.
  }
  \item{degree_centrality_index}{
  The centrality index of the cluster's graph based on within-cluster network degree.
  Computed as the \code{centralization} element of the output from
  \code{\link[igraph:centr_degree]{centr_degree()}}.
  }
  \item{closeness_centrality_index}{
  The centrality index of the cluster's graph based on closeness,
  i.e., distance to other nodes in the cluster.
  Computed using \code{\link[igraph:centralization]{centralization()}}.
  }
  \item{eigen_centrality_index}{
  The centrality index of the cluster's graph based on the eigenvector centrality scores,
  i.e., values of the first eigenvector of the adjacency matrix for the cluster.
  Computed as the \code{centralization} element of the output from
  \code{\link[igraph:centr_eigen]{centr_eigen()}}.
  }
  \item{eigen_centrality_eigenvalue}{
  The eigenvalue corresponding to the first eigenvector of the adjacency matrix
  for the cluster. Computed as the \code{value} element of the output from
  \code{\link[igraph:eigen_centrality]{eigen_centrality()}}.
  }
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
\seealso{
\code{\link[=addClusterStats]{addClusterStats()}}
\code{\link[=addClusterMembership]{addClusterMembership()}}
\code{\link[=labelClusters]{labelClusters()}}
}
\examples{
set.seed(42)
toy_data <- simulateToyData()

net <-
  generateNetworkObjects(
    toy_data, "CloneSeq"
  )

net <- addClusterMembership(net)

net$cluster_data <-
  getClusterStats(
    net$node_data,
    net$adjacency_matrix,
    seq_col = "CloneSeq",
    count_col = "CloneCount"
  )
}
