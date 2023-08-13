\name{plotNetworkGraph}
\alias{plotNetworkGraph}
\title{
  Plot the Graph of an Immune Repertoire Network
}
\description{
  Given the \code{\link[igraph]{igraph}} of an immune repertoire network,
  generates a plot of the network graph according to the user specifications.

  For manual generation of plots, the \code{\link{generateNetworkGraphPlots}}
  function is recommended over  \code{plotNetworkGraph}. The arguments of
  \code{\link{generateNetworkGraphPlots}} behave like those of
  \code{\link{buildRepSeqNetwork}} and support generation of multiple plots.
  \code{\link{generateNetworkGraphPlots}} dispatches appropriate calls to
  \code{plotNetworkGraph}, so most users should have no need to call
  \code{plotNetworkGraph} directly.
}
\usage{
plotNetworkGraph(
    igraph,
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
    outfile = NULL,
    pdf_width = 12,
    pdf_height = 8
)
}
\arguments{
  \item{igraph}{
      An object of class \code{\link[igraph]{igraph}}.
}
  \item{plot_title}{
      A character string containing the plot title. Passed to
      \code{\link[ggplot2]{labs}}.
}
  \item{plot_subtitle}{
      A character string containing the plot subtitle. Passed to
      \code{\link[ggplot2]{labs}}.
}
  \item{color_nodes_by}{
      A vector whose length matches the number of nodes in the network.
      The values are used to encode the color of each node. An argument value of
      \code{NULL} (the default) leaves the nodes uncolored. Passed to the color
      aesthetic mapping of \code{\link[ggraph]{geom_node_point}}.
}
  \item{color_scheme}{
    A character string specifying the color scale used to color the nodes.
    \code{"default"} uses default \code{\link[ggplot2]{ggplot}} colors.
    Other options are one of the viridis color scales (e.g., \code{"plasma"},
    \code{"A"} or other valid inputs to the \code{option} argument of
    \code{\link[ggraph]{scale_color_viridis}}) or (for discrete variables) a palette
    from \code{\link[grDevices]{hcl.pals}} (e.g., \code{"RdYlGn"}). Each of the
    viridis color scales can include the suffix \code{"-1"}
    to reverse its direction (e.g., \code{"plasma-1"} or \code{"A-1"}).
}
  \item{color_legend}{
    A logical scalar specifying whether to display the color legend in the plot.
    The default value of \code{"auto"} shows the color legend if \code{color_nodes_by}
    is a continuous variable or a discrete variable with at most 20 distinct values.
}
  \item{color_title}{
    A character string (or \code{NULL}) specifying the title for the color legend.
    Only applicable if  \code{color_nodes_by} is a vector. If
    \code{color_title = "auto"} (the default), the title for the color legend will
    be the name of the vector provided to \code{color_nodes_by}.
}
  \item{edge_width}{
      A numeric scalar specifying the width of the graph edges in the plot.
      Passed to the \code{width} argument of \code{\link[ggraph]{geom_edge_link0}}.
}
  \item{size_nodes_by}{
    A numeric scalar specifying the size of the nodes, or a numeric vector with
    positive entires that encodes the size of each node (and whose length matches
    the number of nodes in the network). Alternatively, an argument value of
    \code{NULL} uses the default \code{\link[ggraph]{ggraph}} size for all nodes.
    Passed to the size aesthetic mapping of \code{\link[ggraph]{geom_node_point}}.
}
  \item{size_title}{
    A character string (or \code{NULL}) specifying the title for the size legend.
    Only applicable if  \code{size_nodes_by} is a vector. If
    \code{size_title = "auto"} (the default),  the title for the color legend
    will be the name of the vector provided to  \code{size_nodes_by}.
}
  \item{node_size_limits}{
    A numeric vector of length 2, specifying the minimum and maximum node size.
    Only applicable if \code{size_nodes_by} is a vector. If
    \code{node_size_limits = NULL}, the default size scale will be used.
}
  \item{outfile}{
    An optional file path for saving the plot as a pdf. If \code{NULL}
    (the default), no pdf will be saved.
}
  \item{pdf_width}{
    Sets the plot width when writing to pdf. Passed to the \code{width} argument
    of the \code{\link[grDevices]{pdf}} function.
}
  \item{pdf_height}{
    Sets the plot height when writing to pdf. Passed to the \code{height} argument
    of the \code{\link[grDevices]{pdf}} function.
}
}
\value{
  A \code{\link[ggraph]{ggraph}} object.
}
\seealso{
\code{\link{generateNetworkGraphPlots}}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/network_visualization.html}{Network Visualization vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData()

# Generate network for data
net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )

# Plot network graph
net_plot <- plotNetworkGraph(
  net$igraph,
  color_nodes_by =
    net$node_data$SampleID,
  color_title = NULL,
  color_scheme = "viridis",
  size_nodes_by =
    net$node_data$CloneCount,
  size_title = "Clone Count",
  node_size_limits = c(0.5, 3))

print(net_plot)

}