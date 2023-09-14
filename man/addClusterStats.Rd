\name{addClusterStats}
\alias{addClusterStats}
\title{
    Compute Cluster-Level Network Properties
}
\description{
    Given a list of network objects returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}},
    computes cluster-level network properties,
    performing clustering first if needed.
    The list of network objects is returned
    with the cluster properties added as a data frame.
}
\usage{
addClusterStats(
  net,
  cluster_id_name = "cluster_id",
  seq_col = NULL,
  count_col = NULL,
  degree_col = "degree",
  cluster_fun = "fast_greedy",
  overwrite = FALSE,
  verbose = FALSE,
  ...
)

}
\arguments{
  \item{net}{
      A \code{\link{list}} of network objects conforming to the output of
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
      See details.
}
  \item{cluster_id_name}{
      A character string specifying the name of the cluster membership variable
      in \code{net$node_data} that identifies the cluster to which each node belongs.
      If the variable does not exist, it will be added by calling
      \code{\link[=addClusterMembership]{addClusterMembership()}}.
      If the variable does exist, its values will be used
      unless \code{overwrite = TRUE},
      in which case its values will be overwritten and the new values used.
}
  \item{seq_col}{
      Specifies the column(s) of \code{net$node_data} containing
      the receptor sequences upon whose similarity the network is based.
      Accepts a character or numeric vector of length 1 or 2,
      containing either column names or column indices.
      If provided, related cluster-level properties will be computed.
      The default \code{NULL} will use the value contained in
      \code{net$details$seq_col} if it exists and is valid.
}
  \item{count_col}{
      Specifies the column of \code{net$node_data} containing
      a measure of abundance (such as clone count or UMI count).
      Accepts a character string containing the column name
      or a numeric scalar containing the column index.
      If provided,  related cluster-level properties will be computed.
}
  \item{degree_col}{
      Specifies the column of \code{net$node_data} containing
      the network degree of each node.
      Accepts a character string containing the column name.
      If the column does not exist, it will be added.
}
  \item{cluster_fun}{
      A character string specifying the clustering algorithm to use when
      adding or overwriting the cluster membership variable in
      \code{net$node_data} specified by \code{cluster_id_name}.
      Passed to \code{\link[=addClusterMembership]{addClusterMembership()}}.
}
  \item{overwrite}{
      Logical. If \code{TRUE} and \code{net} already contains an element named
      \code{cluster_data}, it will be overwritten. Similarly, if
      \code{overwrite = TRUE} and \code{net$node_data} contains a variable whose
      name matches the value of \code{cluster_id_name}, then its values will be
      overwritten with new cluster membership values (obtained using
      \code{\link[=addClusterMembership]{addClusterMembership()}}
      with the specified value of \code{cluster_fun}),
      and cluster properties will be computed based on the new values.
  }
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
  \item{...}{
      Named optional arguments to the function specified by \code{cluster_fun}.
  }

}
\details{
    The list \code{net} must contain the named elements
    \code{igraph} (of class \code{\link[igraph]{igraph}}),
    \code{adjacency_matrix} (a \code{\link{matrix}} or
    \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} encoding edge connections),
    and \code{node_data} (a \code{\link{data.frame}} containing node metadata),
    all corresponding to the same network. The lists returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} and
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}
    are examples of valid inputs for the \code{net} argument.

    If the network graph has previously been partitioned into clusters using
    \code{\link[=addClusterMembership]{addClusterMembership()}} and the user
    wishes to compute network properties for these clusters, the name of the
    cluster membership variable in \code{net$node_data} should be provided to
    the \code{cluster_id_name} argument.

    If the value of \code{cluster_id_name} is not the name of a variable
    in \code{net$node_data}, then clustering is performed using
    \code{\link[=addClusterMembership]{addClusterMembership()}}
    with the specified value of \code{cluster_fun},
    and the cluster membership values are written to \code{net$node_data} using
    the value of \code{cluster_id_name} as the variable name.
    If \code{overwrite = TRUE}, this is done even if this variable already exists.

}
\value{
  A modified copy of \code{net}, with cluster properties contained in the element
  \code{cluster_data}. This is a \code{\link{data.frame}} containing
  one row for each cluster in the network and the following variables:

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

  If \code{net$node_data} did not previously contain a variable whose name matches
  the value of \code{cluster_id_name}, then this variable will be present
  and will contain values for cluster membership, obtained through a call to
  \code{\link[=addClusterMembership]{addClusterMembership()}}
  using the clustering algorithm specified by \code{cluster_fun}.

  If \code{net$node_data} did previously contain a variable whose name matches
  the value of \code{cluster_id_name} and \code{overwrite = TRUE}, then the
  values of this variable will be overwritten with new values for cluster membership,
  obtained as above based on \code{cluster_fun}.

  If \code{net$node_data} did not previously contain a variable whose name matches
  the value of \code{degree_col}, then this variable will be present
  and will contain values for network degree.

  Additionally, if \code{net} contains a list named \code{details}, then the
  following elements will be added to \code{net$details}, or overwritten if they
  already exist:

  \item{\code{cluster_data_goes_with}}{
      A character string containing the value of \code{cluster_id_name}.
      When \code{net$node_data} contains multiple cluster membership variables
      (e.g., from applying different clustering methods),
      \code{cluster_data_goes_with} allows the user
      to distinguish which of these variables corresponds to
      \code{net$cluster_data}.
  }
  \item{\code{count_col_for_cluster_data}}{
      A character string containing the value of \code{count_col}.
      If \code{net$node_data} contains multiple count variables, this allows
      the user to distinguish which of these variables corresponds to the
      count-related properties in \code{net$cluster_data}, such as
      \code{max_count}.
      If \code{count_col = NULL}, then the value will be \code{NA}.
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
\code{\link[=addClusterMembership]{addClusterMembership()}}
\code{\link[=getClusterStats]{getClusterStats()}}
\code{\link[=labelClusters]{labelClusters()}}
}
\examples{
set.seed(42)
toy_data <- simulateToyData()

net <- generateNetworkObjects(
  toy_data, "CloneSeq"
)

net <- addClusterStats(
  net,
  count_col = "CloneCount"
)

head(net$cluster_data)
net$details

# won't change net since net$cluster_data exists
net <- addClusterStats(
  net,
  count_col = "CloneCount",
  cluster_fun = "leiden",
  verbose = TRUE
)

# overwrites values in net$cluster_data
# and cluster membership values in net$node_data$cluster_id
# with values obtained using "cluster_leiden" algorithm
net <- addClusterStats(
  net,
  count_col = "CloneCount",
  cluster_fun = "leiden",
  overwrite = TRUE
)

net$details

# overwrites existing values in net$cluster_data
# with values obtained using "cluster_louvain" algorithm
# saves cluster membership values to net$node_data$cluster_id_louvain
# (net$node_data$cluster_id retains membership values from "cluster_leiden")
net <- addClusterStats(
  net,
  count_col = "CloneCount",
  cluster_fun = "louvain",
  cluster_id_name = "cluster_id_louvain",
  overwrite = TRUE
)

net$details

# perform clustering using "cluster_fast_greedy" algorithm,
# save cluster membership values to net$node_data$cluster_id_greedy
net <- addClusterMembership(
  net,
  cluster_fun = "fast_greedy",
  cluster_id_name = "cluster_id_greedy"
)

# compute cluster properties for the clusters from previous step
# overwrites values in net$cluster_data
net <- addClusterStats(
  net,
  cluster_id_name = "cluster_id_greedy",
  overwrite = TRUE
)

net$details
}