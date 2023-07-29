\name{addClusterMembership}
\alias{addClusterMembership}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Perform Cluster Analysis on the Nodes of a Network Graph
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given the \code{igraph} and node-level metadata for a network, partitions the graph into densely-connected subgraphs (clusters) using the specified clustering algorithm, and adds a cluster membership ID variable to the node-level metadata.
}
\usage{
addClusterMembership(
  data,
  net,
  fun = cluster_fast_greedy
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
  \item{fun}{
%%     ~~Describe \code{net} here~~
      The function to use for the clustering algorithm. See \code{\link{clustering_algorithms}}.
}
}
\details{
    Cluster analysis involves using a community-finding algorithm to partition the network graph into clusters (densely-connected subgraphs). These clusters represent groups of clones/cells with similar receptor sequences.

    Each cluster is assigned a numeric cluster ID, and the cluster membership ID of each node is recorded. The cluster membership ID is the cluster ID of the cluster to which the node belongs.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A copy of \code{data} with an additional column named \code{cluster_id} that contains the cluster membership ID for the network node corresponding to each row.
}
\seealso{
\code{\link{addClusterLabels}}
\code{\link{getClusterStats}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
