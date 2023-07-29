\name{addNodeNetworkStats}
\alias{addNodeNetworkStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Compute Node-Level Network Properties
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the node-level metadata for the network, such as the data frame named \code{node_data} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{net}{
%%     ~~Describe \code{net} here~~
      The \code{\link[igraph]{igraph}} object for the network, such as the element named \code{igraph} contained in the output of \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}.
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
      Specifies which network properties to compute. Either a vector created using \code{\link{chooseNodeStats}} or \code{\link{exclusiveNodeStats}}, or \code{"all"} to compute all network properties, or \code{"cluster_id_only"} to compute only the \code{cluster_id} property.
}
  \item{cluster_fun}{
%%     ~~Describe \code{net} here~~
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A copy of \code{data} with with an additional column for each network property computed. See \code{\link{chooseNodeStats}} for the network property names, which are used as the column names.
}
\seealso{
\code{\link{chooseNodeStats}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
