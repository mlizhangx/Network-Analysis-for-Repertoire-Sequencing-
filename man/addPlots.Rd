\name{addPlots}
\alias{addPlots}
\alias{generateNetworkGraphPlots}
\title{
  Generate Plots of a Network Graph
}
\description{
  Generates one or more
  \code{\link[ggraph]{ggraph}} plots of the network graph according to the user
  specifications.

  \code{addPlots()} accepts and returns a list of network objects, adding
  the plots to the existing list contents. If the list already contains plots,
  the new plots will be created using the same coordinate layout as the
  existing plots.

  \code{generateNetworkGraphPlots()} accepts the network
  \code{\link[igraph]{igraph}} and node metadata,
  and returns a list containing plots.
}
\usage{
addPlots(
  net,
  print_plots = FALSE,
  plot_title = NULL,
  plot_subtitle = "auto",
  color_nodes_by = NULL,
  color_scheme = "default",
  color_legend = "auto",
  color_title = "auto",
  edge_width = 0.1,
  size_nodes_by = 0.5,
  node_size_limits = NULL,
  size_title = "auto",
  verbose = FALSE
)

generateNetworkGraphPlots(
  igraph,
  data,
  print_plots = FALSE,
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
  layout = NULL,
  verbose = FALSE
)
}

\arguments{

\item{net}{
      A \code{\link{list}} of network objects conforming to the output of
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
      See details.
}

\item{igraph}{
  An \code{\link[igraph]{igraph}} object containing the network graph to be plotted.
}

\item{data}{
  A data frame containing the node metadata for the network, with each row
  corresponding to a node.
}
\item{print_plots}{
  A logical scalar; should plots be printed in the \code{R} plotting window?
}

\item{plot_title}{
      A character string containing the plot title.
}

\item{plot_subtitle}{
      A character string containing the plot subtitle. The default value
      \code{"auto"} generates a subtitle describing the settings used to
      construct the network,
      including the distance type and distance cutoff.
}

\item{color_nodes_by}{
  A vector specifying one or more node metadata variables used to encode the
  color of the nodes. One plot is generated for each entry, with each plot
  coloring the nodes according to the variable in the corresponding entry.
  This argument accepts a character vector where each entry is a
  column name of the node metadata.
  If this argument is \code{NULL}, generates a single plot with uncolored nodes.
}

\item{color_scheme}{
  A character string specifying the color scale to use for all plots, or a
  character vector whose length matches that of \code{color_nodes_by}, with each
  entry specifying the color scale for the corresponding plot.
  \code{"default"} specifies the default
  \code{\link[ggplot2:ggplot]{ggplot()}} color scale.
  Other options are one of the viridis color scales (e.g., \code{"plasma"},
  \code{"A"} or other valid inputs to the \code{option} argument of
  \code{\link[ggraph:scale_color_viridis]{scale_color_viridis()}})
  or (for discrete variables) a palette
  from \code{\link[grDevices:hcl.pals]{hcl.pals()}} (e.g., \code{"RdYlGn"}).
  Each of the viridis color scales can include the suffix \code{"-1"}
  to reverse its direction (e.g., \code{"plasma-1"} or \code{"A-1"}).
}

\item{color_legend}{
    A logical scalar specifying whether to display the color legend in plots.
    The default value of \code{"auto"} shows the color legend if nodes are colored
    according to a continuous variable or
    according to a discrete variable with at most 20 distinct values.
}

\item{color_title}{
  A character string specifying the title of the color legend in all plots, or a
  character vector whose length matches that of \code{color_nodes_by}, with each
  entry specifying the title of the color legend in the corresponding plot.
  Only applicable for plots with colored nodes.
  The value \code{"auto"} uses the corresponding value of \code{color_nodes_by}.
}

\item{edge_width}{
      A numeric scalar specifying the width of the graph edges in the plot.
      Passed to the \code{width} argument of
      \code{\link[ggraph:geom_edge_link0]{geom_edge_link0()}}.
}

\item{size_nodes_by}{
  A numeric scalar specifying the size of the nodes in all plots, or the column
  name of a node metadata variable used to encode the size of the nodes in all
  plots. Alternatively, an argument value of \code{NULL} uses the default
  \code{\link[ggraph]{ggraph}} size for all nodes. Passed to the size aesthetic
  mapping of \code{\link[ggraph:geom_node_point]{geom_node_point()}}.
}

\item{node_size_limits}{
    A numeric vector of length 2, specifying the minimum and maximum node size.
    Only applicable if nodes are sized according to a variable.
    If \code{node_size_limits = NULL}, the default size scale will be used.
}

\item{size_title}{
    A character string (or \code{NULL}) specifying the title for the size legend.
    Only applicable if nodes are sized according to a variable.
    The value \code{"auto"} uses the value of \code{size_nodes_by}.
}
\item{layout}{
  A \code{matrix} specifying the coordinate layout of the network nodes,
  with one row for each node in the network and two columns.
  Each row specifies the x and y coordinates for the corresponding node.
  If \code{NULL}, the layout matrix is created using
  \code{[igraph:layout_components]{layout_components()}}.
  This argument can be used to
  create plots conforming to the same layout as previously-generated plots. It
  can also be used to generate plots with custom layouts.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
}

\details{
    The list \code{net} must contain the named elements
    \code{igraph} (of class \code{\link[igraph]{igraph}}),
    \code{adjacency_matrix} (a \code{\link{matrix}} or
    \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} encoding edge connections),
    and \code{node_data} (a \code{\link{data.frame}} containing node metadata),
    all corresponding to the same network. The lists returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} and
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}
    are examples of valid inputs for the \code{net} argument.

    The arguments \code{color_nodes_by} and \code{size_nodes_by} accept
    the names of variables in the node metadata.
    For \code{addPlots()}, this is the data frame \code{node_data}
    contained in the list provided to the \code{net} argument.
    For \code{generateNetworkGraphPlots()}, this is the data frame provided
    to the \code{data} argument.

    \code{addPlots()} adds the generated plots to the list \code{plots} contained
    in the list of network objects provided to \code{net}.
    The \code{plots} element is created if it does not already exist.
    If plots already exist, the new plots will be generated with the same
    coordinate layout as the existing plots.
    Each plot is named according to the variable used to color the nodes.
    If a plot already exists with the same name as one of the new plots,
    it will be overwritten with the new plot.
    If the \code{plots} list does not already contain an element named
    \code{graph_layout}, it will be added. This element contains the coordinate
    layout for the plots as a two-column matrix.

    When calling \code{generateNetworkGraphPlots()}, if one wishes for the plots
    to be generated with the same coordinate layout as an existing plot, the
    layout matrix for the existing plot must be passed to the \code{layout}
    argument.

    The plots can be printed to a pdf using
    \code{\link[=saveNetworkPlots]{saveNetworkPlots()}}.
}

\value{
  \code{addPlots()} returns a modified copy of \code{net} with the new plots
  contained in the element named \code{plots} (a list), in addition to any
  previously existing plots.

  \code{generateNetworkGraphPlots()} returns a list containing the new plots.

  Each plot is an object of class \code{ggraph}. Within the list of plots,
  each plot is named after the variable used to color the nodes.
  For a plot with uncolored nodes, the name is \code{uniform_color}.

  The list containing the new plots also contains an element named
  \code{graph_layout}. This is
  a matrix specifying the coordinate layout of the nodes in the plots.
  It contains one row for each node in the
  network and two columns. Each row specifies the x and y coordinates for the
  corresponding node. This matrix can be used to generate additional plots with
  the same layout as the plots in the returned list.
}

\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/network_visualization.html}{Network Visualization article on package website}
}

\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}

\seealso{
  \code{\link[=labelNodes]{labelNodes()}}
  \code{\link[=labelClusters]{labelClusters()}}
  \code{\link[=saveNetworkPlots]{saveNetworkPlots()}}
}

\examples{
set.seed(42)
toy_data <- simulateToyData()

net <- buildNet(toy_data, "CloneSeq", node_stats = TRUE)

net <- addPlots(
  net,
  color_nodes_by =
    c("SampleID", "transitivity", "coreness"),
  color_scheme =
    c("Set 2", "mako-1", "plasma-1"),
  color_title =
    c("", "Transitvity", "Coreness"),
  size_nodes_by = "degree",
  node_size_limits = c(0.1, 1.5),
  plot_subtitle = NULL,
  print_plots = TRUE
)

}
