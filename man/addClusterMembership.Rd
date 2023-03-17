\name{addClusterMembership}
\alias{addClusterMembership}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Compute Cluster Membership of Network Nodes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given immune repertoire sequencing data and the
    \code{igraph} object for the network by sequence similarity,
    identifies dense clusters within the network and augments the data
    with a cluster membership variable.
}
\usage{
addClusterMembership(data, net, fun = cluster_fast_greedy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones/cells
      indexed by row.
}
  \item{net}{
%%     ~~Describe \code{net} here~~
      An \code{igraph} object containing the network on the rows of the data.
}
  \item{fun}{
%%     ~~Describe \code{net} here~~
      A function from the \code{igraph} package that takes a graph as its primary argument, implements a community-detection algorithm and returns a community object. See details below.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The clusters are identified by searching for densely-connected subgraphs (with relatively many edge connections among their nodes), also called communities, within the network.

    The \code{fun} argument controls the algorithm used to detect the community structure in the network. Following is a list of all functions in the \code{igraph} package (as of version 1.3.5) that can be used: \code{cluster_fast_greedy}, \code{cluster_edge_betweenness}, \code{cluster_fluid_communities}, \code{cluster_infomap}, \code{cluster_label_prop}, \code{cluster_leading_eigen}, \code{cluster_leiden}, \code{cluster_louvain}, \code{cluster_optimal}, \code{cluster_spinglass},  \code{cluster_walktrap}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A copy of \code{data} with an additional column named \code{cluster_id} that contains the cluster membership for each row.
}
\references{
%% ~put references to the literature/web site here ~
%%    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")

# Add cluster ID
data_w_clusterID <- addClusterMembership(net$node_data, net$igraph)
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
