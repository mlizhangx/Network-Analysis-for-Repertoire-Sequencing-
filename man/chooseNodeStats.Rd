\name{chooseNodeStats}
\alias{chooseNodeStats}
\alias{exclusiveNodeStats}
\alias{node_stat_settings}
\title{
    Specify Node-level Network Properties to Compute
}
\description{
    Create a vector specifying node-level network properties to compute.
    Intended for use with \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
    \code{\link{addNodeNetworkStats}}.

    \code{node_stat_settings()} is a deprecated equivalent of
    \code{chooseNodeStats()}.
}
\usage{
chooseNodeStats(
  degree = TRUE,
  cluster_id = FALSE,
  transitivity = TRUE,
  closeness = FALSE,
  centrality_by_closeness = FALSE,
  eigen_centrality = TRUE,
  centrality_by_eigen = TRUE,
  betweenness = TRUE,
  centrality_by_betweenness = TRUE,
  authority_score = TRUE,
  coreness = TRUE,
  page_rank = TRUE,
  all_stats = FALSE
)

exclusiveNodeStats(
  degree = FALSE,
  cluster_id = FALSE,
  transitivity = FALSE,
  closeness = FALSE,
  centrality_by_closeness = FALSE,
  eigen_centrality = FALSE,
  centrality_by_eigen = FALSE,
  betweenness = FALSE,
  centrality_by_betweenness = FALSE,
  authority_score = FALSE,
  coreness = FALSE,
  page_rank = FALSE
)

}
\arguments{
  \item{degree}{
    Logical. Whether to compute network degree.
}
  \item{cluster_id}{
    Logical. Whether to perform cluster analysis and record the cluster
    membership of each node.
    See \code{\link[=addClusterMembership]{addClusterMembership()}}.
}
  \item{transitivity}{
    Logical. Whether to compute node-level network transitivity using
    \code{\link[igraph:transitivity]{transitivity()}}
    with \code{type = "local"}. The local
    transitivity of a node is the the number of triangles connected to the node
    relative to the number of triples centered on that node.
}
  \item{closeness}{
    Logical. Whether to compute network closeness using
    \code{\link[igraph:closeness]{closeness()}}.
}
  \item{centrality_by_closeness}{
    Logical. Whether to compute network centrality by closeness. The values are
    the entries of the \code{res} element of the list returned by
    \code{\link[igraph:centr_clo]{centr_clo()}}.
}
  \item{eigen_centrality}{
    Logical. Whether to compute the eigenvector centrality scores of node network
    positions. The scores are the entries of the \code{vector} element of the list
    returned by \code{\link[igraph:eigen_centrality]{eigen_centrality()}}
    with \code{weights = NA}.
    The centrality scores correspond to the values of the first eigenvector of
    the adjacency matrix for the cluster graph.
}
  \item{centrality_by_eigen}{
    Logical. Whether to compute node-level network centrality scores based on
    eigenvector centrality scores. The scores are the entries of the \code{vector}
    element of the list returned by \code{\link[igraph:centr_eigen]{centr_eigen()}}.
}
  \item{betweenness}{
    Logical. Whether to compute network betweenness using
    \code{\link[igraph:betweenness]{betweenness()}}.
}
  \item{centrality_by_betweenness}{
    Logical. Whether to compute network centrality scores by betweenness. The
    scores are the entires of the \code{res} element of the list returned by
    \code{\link[igraph:centr_betw]{centr_betw()}}.
}
  \item{authority_score}{
    Logical. Whether to compute the authority score using
    \code{\link[igraph:authority_score]{authority_score()}}.
}
  \item{coreness}{
    Logical. Whether to compute network coreness using
    \code{\link[igraph:coreness]{coreness()}}.
}
  \item{page_rank}{
    Logical. Whether to compute page rank. The page rank values are the entries of
    the \code{vector} element of the list returned by
    \code{\link[igraph:page_rank]{page_rank()}}.
}
  \item{all_stats}{
    Logical. If \code{TRUE}, all other argument values are overridden and set to
    \code{TRUE}.
}
}
\details{
    These functions return a vector that can be passed to the \code{stats_to_include}
    argument of \code{\link[=addNodeStats]{addNodeStats()}}
    (or \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}, if
    \code{node_stats = TRUE})
    in order to specify which node-level network properties to compute.

    \code{chooseNodeStats} and \code{exclusiveNodeStats} each have default
    argument values suited to a different use case,
    in order to reduce the number of argument values that must be set manually.

    \code{chooseNodeStats} has most arguments \code{TRUE} by default.
    It is best suited for including a majority of the available properties.
    It can be called with \code{all_stats = TRUE} to set all values to \code{TRUE}.

    \code{exclusiveNodeStats} has all of its arguments set to \code{FALSE} by
    default. It is best suited for including only a few properties.
}
\value{
    A named logical vector with one entry for each of the function's arguments
    (except for \code{all_stats}).
    Each entry has the same name as the corresponding argument, and its value
    matches the argument's value.
}
\seealso{
\code{\link[=addNodeStats]{addNodeStats()}}
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
