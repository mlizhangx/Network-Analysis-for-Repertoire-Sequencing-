\name{node_stat_settings}
\alias{node_stat_settings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Specify Node Network Statistics to Compute
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Create a list, for use by \code{addNodeNetworkStats()}, specifying the node-level
    network stats to be computed.
}
\usage{
node_stat_settings(degree = TRUE,
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
                   all_stats = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{degree}{
%%     ~~Describe \code{degree} here~~
    Logical. Whether to compute network degree.
}
  \item{cluster_id}{
%%     ~~Describe \code{cluster_id} here~~
    Logical. Whether to perform clustering, which is done using \code{igraph::cluster_fast_greedy()}, in order to determine the cluster membership of each node.
}
  \item{transitivity}{
%%     ~~Describe \code{transitivity} here~~
    Logical. Whether to compute node-level network transitivity, which is done using \code{igraph::transitivity()}
    with \code{type = "local"} The local transitivity of a node is the the number of triangles connected to the node relative to the number of triples centered on that node.
}
  \item{closeness}{
%%     ~~Describe \code{closeness} here~~
    Logical. Whether to compute network closeness, which is done using \code{igraph::closeness()}.
}
  \item{centrality_by_closeness}{
%%     ~~Describe \code{centrality_by_closeness} here~~
    Logical. Whether to compute network centrality by closeness, which is done using \code{igraph::centr_clo()$res}.
}
  \item{eigen_centrality}{
%%     ~~Describe \code{eigen_centrality} here~~
    Logical. Whether to compute the eigenvector centrality scores of node network positions, which is done using \code{igraph::eigen_centrality()$vector} with \code{weights = NA}. The centrality scores correspond to the values of the first eigenvector of the adjacency matrix for the cluster graph.
}
  \item{centrality_by_eigen}{
%%     ~~Describe \code{centrality_by_eigen} here~~
    Logical. Whether to compute node-level network centrality scores based on eigenvector centrality scores, which is done using \code{igraph::centr_eigen()$vector}.
}
  \item{betweenness}{
%%     ~~Describe \code{betweenness} here~~
    Logical. Whether to compute network betweenness, which is done using \code{igraph::betweenness()}.
}
  \item{centrality_by_betweenness}{
%%     ~~Describe \code{centrality_by_betweenness} here~~
    Logical. Whether to compute network centrality by betweenness, which is done using \code{igraph::centr_betw()$res}.
}
  \item{authority_score}{
%%     ~~Describe \code{authority_score} here~~
    Logical. Whether to compute the authority score, which is done using \code{igraph::authority_score()}.
}
  \item{coreness}{
%%     ~~Describe \code{coreness} here~~
    Logical. Whether to compute network coreness, which is done using \code{igraph::coreness()}.
}
  \item{page_rank}{
%%     ~~Describe \code{page_rank} here~~
    Logical. Whether to compute page rank, which is done using \code{igraph::page_rank()$vector}.
}
  \item{all_stats}{
%%     ~~Describe \code{all_stats} here~~
    Logical. If \code{TRUE}, then \code{addNodeNetworkStats()} will disregard the values for the other settings and instead treat all of them as \code{TRUE}.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    An object of type \code{"list"}, for use with \code{addNodeNetworkStats()} as input to the \code{stats_to_include} argument.
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
\dontrun{
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")

# Add custom network statistics
data_w_stats <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include =
      node_stat_settings(
          cluster_id = TRUE,
          closeness = TRUE,
          centrality_by_closeness = TRUE,
          betweenness = FALSE,
          centrality_by_betweenness = FALSE,
          authority_score = FALSE,
          page_rank = FALSE
      )
  )

# Add all network statistics
data_w_all_stats <-
  addNodeNetworkStats(
    net$node_data, net$igraph,
    stats_to_include =
      node_stat_settings(all_stats = TRUE)
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
