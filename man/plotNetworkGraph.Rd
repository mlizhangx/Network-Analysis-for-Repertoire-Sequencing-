\name{plotNetworkGraph}
\alias{plotNetworkGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Plot the Graph of an Immune Repertoire Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given an \code{\link[igraph]{igraph}} network object, generates a \code{\link[ggraph]{ggraph}} plot of the network graph according to the user specifications.

  For manual generation of plots, the \code{\link{generateNetworkGraphPlots}} function is recommended over  \code{plotNetworkGraph}. The arguments of \code{\link{generateNetworkGraphPlots}} behave like those of \code{\link{buildRepSeqNetwork}} and support generation of multiple plots. \code{\link{generateNetworkGraphPlots}} dispatches appropriate calls to \code{plotNetworkGraph}, so most users should have no need to call \code{plotNetworkGraph} directly.
}
\usage{
plotNetworkGraph(igraph,
                 plot_title = NULL,
                 plot_subtitle = NULL,
                 color_nodes_by = NULL,
                 color_scheme = "default",
                 color_legend = "auto",
                 color_title = "auto",
                 edge_width = 0.1,
                 size_nodes_by = 0.5,
                 node_size_limits = NULL,
                 size_title = "auto",
                 outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{igraph}{
%%     ~~Describe \code{network} here~~
      An \code{\link[igraph]{igraph}} object containing the network graph to be plotted.
}
  \item{plot_title}{
%%     ~~Describe \code{title} here~~
      A character string or \code{NULL}. The plot title, to be passed to \code{\link[ggplot2]{labs}}.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
      A character string or \code{NULL}. The plot subtitle, to be passed to \code{\link[ggplot2]{labs}}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
      A vector whose length matches the number of nodes in the network. The values are used to encode the color of each node. Alternatively, an argument value of \code{NULL} (the default) leaves the nodes uncolored. Passed to the color aesthetic mapping of \code{\link[ggraph]{geom_node_point}}.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    A character string specifying the color scale used to color the nodes. \code{"default"} uses default \code{\link[ggplot2]{ggplot}} colors. Other options are one of the \code{\link[viridisLite]{viridis}} color scales (e.g. \code{"plasma"} or \code{"A"}) or (for discrete variables) a palette from \code{\link[grDevices]{hcl.pals}} (e.g. \code{"RdYlGn"}).
}
  \item{color_legend}{
%%     ~~Describe \code{color_title} here~~
    A logical scalar specifying whether to display the color legend on the plot. The default value of \code{"auto"} shows the color legend if \code{color_nodes_by} is a continuous variable or a discrete variable with at most 20 distinct values.
}
  \item{color_title}{
%%     ~~Describe \code{color_title} here~~
    A character string (or \code{NULL}) specifying the title for the color legend. Only applicable if  \code{color_nodes_by} is a vector. If \code{color_title = "auto"} (the default), the color legend title will use the name of the \code{color_nodes_by} vector.
}
  \item{edge_width}{
%%     ~~Describe \code{edge_width} here~~
      A numeric scalar specifying the width of the graph edges in the plot. Passed to the \code{width} argument of \code{\link[ggraph]{geom_edge_link0}}.
}
  \item{size_nodes_by}{
%%     ~~Describe \code{size_nodes_by} here~~
    A numeric scalar specifying the size of the nodes, or a numeric vector with positive entires that encodes the size of each node (and whose length matches the number of nodes in the network). Alternatively, an argument value of \code{NULL} uses the default \code{\link[ggraph]{ggraph}} size for all nodes. Passed to the size aesthetic mapping of \code{\link[ggraph]{geom_node_point}}.
}
  \item{size_title}{
%%     ~~Describe \code{size_title} here~~
    A character string (or \code{NULL}) specifying the title for the size legend. Only applicable if  \code{size_nodes_by} is a vector. If \code{size_title = "auto"} (the default), the size legend title will use the name of the \code{size_nodes_by} vector.
}
  \item{node_size_limits}{
%%     ~~Describe \code{node_size_limits} here~~
    A numeric vector of length 2 specifying the minimum and maximum node size. Only applicable if \code{size_nodes_by} is a vector. If \code{node_size_limits = NULL}, the default size scale will be used.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~
    An optional file path for saving the plot as a pdf. If \code{NULL} (the default), no pdf will be saved.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{\link[ggraph]{ggraph}} object.
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
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")

# Plot network graph
net_plot <- plotNetworkGraph(
  net$igraph,
  color_nodes_by = net$node_data$SampleID,
  color_scheme = "viridis",
  size_nodes_by = net$node_data$CloneCount,
  node_size_limits = c(0.5, 3))

print(net_plot)
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
