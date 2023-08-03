\name{chooseNodeStats}
\alias{chooseNodeStats}
\alias{exclusiveNodeStats}
\alias{node_stat_settings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Specify Node-level Network Properties to Compute
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Create a vector specifying node-level network properties to compute. Intended for use with \code{\link{buildRepSeqNetwork}} or \code{\link{addNodeNetworkStats}}.
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{degree}{
%%     ~~Describe \code{degree} here~~
    Logical. Whether to compute network degree.
}
  \item{cluster_id}{
%%     ~~Describe \code{cluster_id} here~~
    Logical. Whether to perform cluster analysis and record the cluster membership ID of each node. See \code{\link{addClusterMembership}}.
}
  \item{transitivity}{
%%     ~~Describe \code{transitivity} here~~
    Logical. Whether to compute node-level network transitivity using \code{\link[igraph]{transitivity}} with \code{type = "local"}. The local transitivity of a node is the the number of triangles connected to the node relative to the number of triples centered on that node.
}
  \item{closeness}{
%%     ~~Describe \code{closeness} here~~
    Logical. Whether to compute network closeness using \code{\link[igraph]{closeness}}.
}
  \item{centrality_by_closeness}{
%%     ~~Describe \code{centrality_by_closeness} here~~
    Logical. Whether to compute network centrality by closeness. The values are the entries of the \code{res} element of the list returned by \code{\link[igraph]{centr_clo}}.
}
  \item{eigen_centrality}{
%%     ~~Describe \code{eigen_centrality} here~~
    Logical. Whether to compute the eigenvector centrality scores of node network positions. The scores are the entries of the \code{vector} element of the list returned by \code{\link[igraph]{eigen_centrality}} with \code{weights = NA}. The centrality scores correspond to the values of the first eigenvector of the adjacency matrix for the cluster graph.
}
  \item{centrality_by_eigen}{
%%     ~~Describe \code{centrality_by_eigen} here~~
    Logical. Whether to compute node-level network centrality scores based on eigenvector centrality scores. The scores are the entries of the \code{vector} element of the list returned by \code{\link[igraph]{centr_eigen}}.
}
  \item{betweenness}{
%%     ~~Describe \code{betweenness} here~~
    Logical. Whether to compute network betweenness using \code{\link[igraph]{betweenness}}.
}
  \item{centrality_by_betweenness}{
%%     ~~Describe \code{centrality_by_betweenness} here~~
    Logical. Whether to compute network centrality scores by betweenness. The scores are the entires of the \code{res} element of the list returned by \code{\link[igraph]{centr_betw}}.
}
  \item{authority_score}{
%%     ~~Describe \code{authority_score} here~~
    Logical. Whether to compute the authority score using \code{\link[igraph]{authority_score}}.
}
  \item{coreness}{
%%     ~~Describe \code{coreness} here~~
    Logical. Whether to compute network coreness using \code{\link[igraph]{coreness}}.
}
  \item{page_rank}{
%%     ~~Describe \code{page_rank} here~~
    Logical. Whether to compute page rank. The page rank values are the entries of the \code{vector} element of the list returned by \code{\link[igraph]{page_rank}}.
}
  \item{all_stats}{
%%     ~~Describe \code{all_stats} here~~
    Logical. If \code{TRUE}, then all node-level network properties will be computed regardless of the values used for the other arguments.
}
}
\details{
    These functions return a vector that can be passed to the \code{stats_to_include} argument of \code{\link{addNodeNetworkStats}} or \code{\link{buildRepSeqNetwork}} in order to specify which node-level network properties to compute. Note that \code{\link{buildRepSeqNetwork}} must be called with \code{node_stats = TRUE} for this to be relevant.

    The main difference between \code{chooseNodeStats} and \code{exclusiveNodeStats} lies in their default argument values. The defaults for each are suited to different use cases in order to reduce the number of argument values that must be set manually.

    \code{chooseNodeStats} has most of its arguments set to \code{TRUE} by default. Thus it is most suitable  when the user wishes to include a majority of the available properties. It can also be called with \code{all_stats = TRUE} as a convenient way of computing all available properties.

    \code{exclusiveNodeStats} has all of its arguments set to \code{FALSE} by default. It is therefore most suitable when the user wishes to include only a few properties.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A named logical vector with one entry for each of the function's arguments. Each entry has the same name as the corresponding argument, and its value matches the argument's value.
}
\seealso{
\code{\link{addNodeNetworkStats}}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. https://doi.org/10.3389/fimmu.2023.1181825

\url{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825}

\url{https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-}
}
\author{
%%  ~~who you are~~
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
    net$node_data, net$igraph,
    stats_to_include =
      chooseNodeStats()
  )
# Equivalent form for the above
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
    stats_to_include =
      chooseNodeStats(
        all_stats = TRUE
      )
  )
# Equivalent form for the above
data_w_all_stats <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include = "all"
  )

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
