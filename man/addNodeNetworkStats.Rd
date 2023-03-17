\name{addNodeNetworkStats}
\alias{addNodeNetworkStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Compute Node-Level Network Statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given immune repertoire sequence data and the corresponding repertoire network,
  compute the specified node-level network statistics, which are appended as
  additional columns to the rep-seq data before it is returned.

}
\usage{
addNodeNetworkStats(data, net,
                    stats_to_include = node_stat_settings(),
                    cluster_fun = cluster_fast_greedy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones
      indexed by row.
}
  \item{net}{
%%     ~~Describe \code{net} here~~
      An \code{igraph} object containing the network corresponding to the
      clones in \code{data}.
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
      A list generated using \code{node_stat_settings()} that specifies the network
      statistics to compute. Alternatively, \code{stats_to_include = "all"} will
      compute all node-level network stats, and \code{stats_to_include = "cluster_id_only"}
      will compute cluster membership only.
}
  \item{cluster_fun}{
%%     ~~Describe \code{net} here~~
      Passed to \code{addClusterMembership}. Controls the algorithm used for cluster identification.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A data frame with the same rows and columns as \code{data}, with additional columns
  containing the network statistics for the clone in each row.
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

# Add default network statistics
data <- addNodeNetworkStats(net$node_data, net$igraph)
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
