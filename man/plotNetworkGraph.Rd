\name{plotNetworkGraph}
\alias{plotNetworkGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Plot the Graph of an Immune Repertoire Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  A wrapper function for customized calling of functions from the \code{ggraph} package.
  Given an \code{igraph} network object, generates a \code{ggraph} according to the user
  specifications.
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
      An \code{igraph} object containing the network graph to be plotted.
}
  \item{plot_title}{
%%     ~~Describe \code{title} here~~
      A character string or \code{NULL}. The plot title, to be passed to \code{ggplot2::labs}.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
      A character string or \code{NULL}. The plot subtitle, to be passed to \code{ggplot2::labs}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
      Passed to the color aesthetic mapping of \code{ggraph::geom_node_point}, used to encode the color of each node. Either a vector whose length matches the number of nodes in \code{network}, or \code{NULL}. If the vector is numeric and contains non-integer values, a continuous color scale will be used; otherwise, a discrete color scale will be used.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    A character string specifying the color palette used to color the nodes: either "default" for default \code{ggplot2} colors, a viridis color map option (e.g., \code{"viridis"}, \code{"plasma"}, etc., or one the corresponding strings "A" through "H"; see \code{?viridis} for details), or a palette from \code{grDevices::hcl.pals()} (these can only be used with discrete color scales).
}
  \item{color_legend}{
%%     ~~Describe \code{color_title} here~~
    Accepts a logical scalar specifying whether to display the color legend on the plot. The default value of \code{"auto"} shows the color legend if \code{color_nodes_by} is a continuous variable or a discrete variable with at most 20 distinct values.
}
  \item{color_title}{
%%     ~~Describe \code{color_title} here~~
    A character string or \code{NULL} specifying the title for the color legend. If "auto" and a vector was supplied for \code{color_nodes_by}, the color legend title will attempt to use the name of the supplied vector.
}
  \item{edge_width}{
%%     ~~Describe \code{edge_width} here~~
      A numeric scalar specifying the width of the graph edges in the plot. Passed to the \code{width} argument of \code{ggraph::geom_edge_link0}.
}
  \item{size_nodes_by}{
%%     ~~Describe \code{size_nodes_by} here~~
    Passed to the size aesthetic mapping of \code{ggraph::geom_node_point}, used to encode the size of each node. Either a numeric vector with positive entires whose length matches the number of nodes in \code{network}, a positive-valued number (denoting a fixed constant size for each node), or \code{NULL}.
}
  \item{size_title}{
%%     ~~Describe \code{size_title} here~~
    A character string or \code{NULL} specifying the title for the size legend. The size legend is only displayed if a vector was supplied for \code{size_nodes_by}, in which case \code{size_title = "auto"} will attempt to use the name of the supplied vector.
}
  \item{node_size_limits}{
%%     ~~Describe \code{node_size_limits} here~~
    A numeric vector of length 2 with positive entries, or \code{NULL}. If a vector of values is provided, the node sizes will be rescaled using the first value for the minimum node size and the second value for the maximum node sizes. Only applicable if a vector is supplied for \code{size_nodes_by}.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~
    A character string containing a valid file name ending in \code{".pdf"}, or \code{NULL}. If a file name is supplied, the plot will be saved to the file.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{ggraph} object.
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
# Generate some data
set.seed(42)
sample_size <- 200
group_labels <- rep("Control", times = sample_size)
assign_group <- sample(c(TRUE, FALSE), size = sample_size, replace = TRUE)
base_seq_length <- 7
characters <- c("G", "A", "T", "C")
char_probs <- c(2/3, 2/9, 1/18, 1/18)
characters_sample <- sample(characters,
                            size = base_seq_length * sample_size,
                            replace = TRUE, prob = char_probs)
characters_sample <- matrix(characters_sample, ncol = sample_size)
clone_seqs <- apply(characters_sample, MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
clone_seqs_append <- rep("AATC", times = sample_size)
append_latent_prob <- runif(sample_size, min = 0, max = 1)
for (i in 1:sample_size) {
  case_group <- assign_group[[i]]
  if (case_group) group_labels[[i]] <- "Case"
  if ((case_group & append_latent_prob[[i]] > 0.9) |
      (!case_group & append_latent_prob[[i]] > 0.5)) {
    clone_seqs_append[[i]] <- "AATCGGGG"
  } else if ((case_group & append_latent_prob[[i]] > 0.8) |
             (!case_group & append_latent_prob[[i]] > 0.3)) {
    clone_seqs_append[[i]] <- "AATCGGT"
  } else if ((case_group & append_latent_prob[[i]] > 0.7) |
             (!case_group & append_latent_prob[[i]] > 0.2)) {
    clone_seqs_append[[i]] <- "AATCGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.5) |
             (!case_group & append_latent_prob[[i]] > 0.1)) {
    clone_seqs_append[[i]] <- "AATTGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.25) |
             (!case_group & append_latent_prob[[i]] > 0.05)) {
    clone_seqs_append[[i]] <- "AATTG"
  }
  clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
                            collapse = "")
}
counts <- rbinom(sample_size, size = 300, prob = 0.1)
frequencies <- counts/sum(counts)
data <- data.frame("clone_seq" = clone_seqs,
                   "count" = counts,
                   "frequency" = frequencies,
                   "group" = group_labels)

# Generate network for data
net <- generateNetworkObjects(data, "clone_seq")

# Plot network graph
net_plot <- plotNetworkGraph(
  net$igraph,
  color_nodes_by = net$node_data$group,
  color_scheme = "viridis",
  size_nodes_by = net$node_data$count,
  node_size_limits = c(0.5, 5))
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
