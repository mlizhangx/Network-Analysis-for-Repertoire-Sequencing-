\name{addNodeStats}
\alias{addNodeStats}
\title{
  Compute Node-Level Network Properties
}
\description{
  Given a list of network objects returned by
  \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
  \code{\link[=generateNetworkObjects]{generateNetworkObjects()}},
  computes a specified set of network properties for the network nodes.
  The list of network objects is returned
  with each property added as a variable to the node metadata.
}
\usage{
addNodeStats(
  net,
  stats_to_include = chooseNodeStats(),
  cluster_fun = "fast_greedy",
  cluster_id_name = "cluster_id",
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
  \item{stats_to_include}{
      Specifies which network properties to compute.
      Accepts a vector created using
      \code{\link[=chooseNodeStats]{chooseNodeStats()}} or
      \code{\link[=exclusiveNodeStats]{exclusiveNodeStats()}},
      or the character string \code{"all"} to compute all network properties.
}
  \item{cluster_fun}{
      A character string specifying the clustering algorithm to use when
      computing cluster membership.
      Applicable only when \code{stats_to_include = "all"} or
      \code{stats_to_include["cluster_id"]} is \code{TRUE}.
      Passed to \code{\link[=addClusterMembership]{addClusterMembership()}}.
}
  \item{cluster_id_name}{
      A character string specifying the name of the cluster membership variable
      to be added to the node metadata.
      Applicable only when \code{stats_to_include = "all"} or
      \code{stats_to_include["cluster_id"]} is \code{TRUE}.
      Passed to \code{\link[=addClusterMembership]{addClusterMembership()}}.
}
  \item{overwrite}{
      Logical. If \code{TRUE} and \code{net$node_data} contains a variable whose
      name matches the value of \code{cluster_id_name}, then its values will be
      overwritten with new cluster membership values (obtained using
      \code{\link[=addClusterMembership]{addClusterMembership()}}, to which
      the values of \code{cluster_fun}, \code{overwrite}).
      Applicable only when \code{stats_to_include = "all"} or
      \code{stats_to_include["cluster_id"]} is \code{TRUE}.
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
  Node-level network properties are properties that pertain to each individual
  node in the network graph.

  Some are local properties, meaning that their value for a given node depends only
  on a subset of the nodes in the network. One example is the network degree of a
  given node, which represents the number of other nodes that are directly joined
  to the given node by an edge connection.

  Other properties are global properties, meaning that their value for a given node
  depends on all of the nodes in the network. An example is the authority score of
  a node, which is computed using the entire graph adjacency matrix (if we denote
  this matrix by \eqn{A}, then the principal eigenvector of \eqn{A^T A} represents
  the authority scores of the network nodes).

  See \code{\link[=chooseNodeStats]{chooseNodeStats()}}
  for a list of the available node-level network properties.

  The list \code{net} must contain the named elements
  \code{igraph} (of class \code{\link[igraph]{igraph}}),
  \code{adjacency_matrix} (a \code{\link{matrix}} or
  \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} encoding edge connections),
  and \code{node_data} (a \code{\link{data.frame}} containing node metadata),
  all corresponding to the same network. The lists returned by
  \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} and
  \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}
  are examples of valid inputs for the \code{net} argument.

}
\value{
  A modified copy of \code{net},
  with \code{net$node_data} containing an additional column
  for each new network property computed.
  See \code{\link[=chooseNodeStats]{chooseNodeStats()}}
  for the network property names,
  which are used as the column names,
  except for the cluster membership variable,
  whose name is the value of \code{cluster_id_name}.
}
\seealso{
\code{\link[=chooseNodeStats]{chooseNodeStats()}}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology},
vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
set.seed(42)
toy_data <- simulateToyData()

net <- generateNetworkObjects(
  toy_data, "CloneSeq"
)

# Add default set of node properties
net <- addNodeStats(net)

# Modify default set of node properties
net <- addNodeStats(
  net,
  stats_to_include =
    chooseNodeStats(
      closeness = TRUE,
      page_rank = FALSE
    )
)

# Add only the spepcified node properties
net <- addNodeStats(
  net,
  stats_to_include =
    exclusiveNodeStats(
      degree = TRUE,
      transitivity = TRUE
    )
)

# Add all node-level network properties
net <- addNodeStats(
  net,
  stats_to_include = "all"
)


}