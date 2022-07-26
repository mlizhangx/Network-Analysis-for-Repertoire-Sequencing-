\name{plotNetworkGraph}
\alias{plotNetworkGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Plot the Graph of an Immune Repertoire Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  A wrapper function for customized calling of functions from the `ggraph` package.
  Given an `igraph` network object, generates a `ggraph` according to the user
  specifications.
}
\usage{
plotNetworkGraph(network,
                 edge_width = 0.3,
                 title = NULL,
                 subtitle = NULL,
                 color_nodes_by = NULL,
                 color_scheme = "default",
                 show_color_legend = TRUE,
                 color_legend_title = NULL,
                 size_nodes_by = NULL,
                 node_size_limits = NULL,
                 size_legend_title = NULL,
                 outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{
%%     ~~Describe \code{network} here~~
      An `igraph` object containing the network graph to be plotted.
}
  \item{edge_width}{
%%     ~~Describe \code{edge_width} here~~
      A numeric scalar specifying the width of the graph edges in the plot. Passed to the `width` argument of `ggraph::geom_edge_link0`.
}
  \item{title}{
%%     ~~Describe \code{title} here~~
      A character string or `NULL`. The plot title, to be passed to `ggplot2::labs`.
}
  \item{subtitle}{
%%     ~~Describe \code{subtitle} here~~
      A character string or `NULL`. The plot subtitle, to be passed to `ggplot2::labs`.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
      Passed to the color aesthetic mapping of `ggraph::geom_node_point`, used to encode the color of each node. Either a vector whose length matches the number of nodes in `network`, or `NULL`. If the vector is numeric and contains non-integer values, a continuous color scale will be used; otherwise, a discrete color scale will be used.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    a fixed color specified by name (e.g., "blue") or hexadecimal RGB triplet (e.g., "#6600CC"),
}
  \item{show_color_legend}{
%%     ~~Describe \code{color_legend_title} here~~
}
  \item{color_legend_title}{
%%     ~~Describe \code{color_legend_title} here~~
}
  \item{size_nodes_by}{
%%     ~~Describe \code{size_nodes_by} here~~
}
  \item{size_legend_title}{
%%     ~~Describe \code{size_legend_title} here~~
}
  \item{node_size_limits}{
%%     ~~Describe \code{node_size_limits} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (network, edge_width = 0.3, title = NULL, subtitle = NULL,
    color_nodes_by = NULL, size_nodes_by = NULL, color_legend_title = NULL,
    size_legend_title = NULL, color_scheme = "default", node_size_limits = NULL,
    outfile = NULL)
{
    set.seed(9999)
    layout <- igraph::layout_components(network)
    graph_plot <- ggraph::ggraph(network, layout = layout) +
        ggraph::geom_edge_link0(width = edge_width, colour = "grey")
    if (!is.null(color_nodes_by)) {
        if (!is.null(size_nodes_by)) {
            if (is.numeric(size_nodes_by) & length(size_nodes_by) ==
                1) {
                graph_plot <- graph_plot + ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by),
                  size = size_nodes_by)
            }
            else if (length(size_nodes_by) > 1) {
                graph_plot <- graph_plot + ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by,
                  size = size_nodes_by))
                if (!is.null(node_size_limits)) {
                  graph_plot <- graph_plot + ggplot2::scale_size(range = node_size_limits)
                }
            }
        }
        else {
            graph_plot <- graph_plot + ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by))
        }
    }
    else {
        if (!is.null(size_nodes_by)) {
            if (is.numeric(size_nodes_by) & length(size_nodes_by) ==
                1) {
                graph_plot <- graph_plot + ggraph::geom_node_point(size = size_nodes_by)
            }
            else if (length(size_nodes_by) > 1) {
                graph_plot <- graph_plot + ggraph::geom_node_point(ggplot2::aes(size = size_nodes_by))
                if (!is.null(node_size_limits)) {
                  graph_plot <- graph_plot + ggplot2::scale_size(range = node_size_limits)
                }
            }
        }
        else {
            graph_plot <- graph_plot + ggraph::geom_node_point()
        }
    }
    graph_plot <- graph_plot + ggraph::theme_graph(base_family = "sans") +
        ggplot2::labs(title = title, subtitle = subtitle) + ggplot2::guides(color = ggplot2::guide_legend(title = color_legend_title),
        size = ggplot2::guide_legend(title = size_legend_title))
    color_type <- ggplot2::scale_type(color_nodes_by)[[1]]
    if (color_scheme != "default") {
        if (color_type == "continuous") {
            if (color_scheme \%in\% c("A", "B", "C", "D", "E",
                "F", "G", "H", "magma", "inferno", "plasma",
                "viridis", "cividis", "rocket", "mako", "turbo")) {
                graph_plot <- graph_plot + ggraph::scale_color_viridis(option = color_scheme)
            }
            else {
                warning("'color_nodes_by' is continuous; 'color_scheme' must be 'default' or a viridis color map option (see `?viridis`); using default color scheme instead")
            }
        }
        else {
            if (color_scheme \%in\% c("A", "B", "C", "D", "E",
                "F", "G", "H", "magma", "inferno", "plasma",
                "viridis", "cividis", "rocket", "mako", "turbo")) {
                graph_plot <- graph_plot + ggraph::scale_color_viridis(option = color_scheme,
                  discrete = TRUE)
            }
            else if (color_scheme \%in\% grDevices::hcl.pals()) {
                graph_plot <- graph_plot + ggplot2::scale_color_manual(values = grDevices::hcl.colors(n = length(color_nodes_by),
                  palette = color_scheme))
            }
            else {
                warning("'color_scheme' must be 'default' or one of the values contained in `grDevices::hcl.pals()`; using default color scheme instead")
            }
        }
    }
    if (!is.null(outfile)) {
        grDevices::pdf(file = outfile, width = 12, height = 8)
        print(graph_plot)
        grDevices::dev.off()
        cat(paste0("Plot of network graph saved to file:\n  ",
            outfile, "\n"))
    }
    return(graph_plot)
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
