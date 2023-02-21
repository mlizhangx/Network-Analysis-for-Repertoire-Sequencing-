\name{addClusterLabels}
\alias{addClusterLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Label Notable Clusters in a Graph Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given a graph plot and clustering analysis output, ranks clusters according to
  specified property (node size by default) and labels those above a certain rank with
  their cluster IDs, returning the labeled plot.
}
\usage{
addClusterLabels(
  plot, net,
  top_n_clusters = 20,
  criterion = "node_count",
  size = 5, color = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plot}{
%%     ~~Describe \code{plot} here~~
  A \code{ggraph} object.
}
  \item{net}{
%%     ~~Describe \code{net} here~~
  A list of network objects containing list items \code{node_data} and \code{cluster_data}, such as the output of \code{buildRepSeq} when called with \code{cluster_stats = TRUE}.
}
  \item{top_n_clusters}{
%%     ~~Describe \code{top_n_clusters} here~~
  A positive integer specifying the number of clusters to label. Those with the highest rank according to the \code{criterion} argument will be labeled.
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
  The name of a numeric cluster property by which to rank the clusters. See \code{\link{getClusterStats}} for a list of properties. The default is \code{"node_size"}.
}
  \item{size}{
%%     ~~Describe \code{size} here~~
  A positive value specifying the size of the text labels.
}
  \item{color}{
%%     ~~Describe \code{color} here~~
  The color of the text labels.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{ggraph} object; the input plot, annotated with cluster ID labels according to the user specifications.
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

# Generate network, but don't print the graph plot
network <- buildRepSeqNetwork(
  toy_data, seq_col = "CloneSeq",
  node_stats = TRUE, cluster_stats = TRUE,
  color_nodes_by = "transitivity", color_scheme = "plasma-1",
  size_nodes_by = "degree", node_size_limits = c(0.5, 1.5),
  print_plots = FALSE, output_dir = NULL)

# Add labels to the two largest clusters and print the plot
addClusterLabels(plot = network$plots$transitivity,
                 net = network,
                 top_n_clusters = 2,
                 criterion = "node_count" # (the default)
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
