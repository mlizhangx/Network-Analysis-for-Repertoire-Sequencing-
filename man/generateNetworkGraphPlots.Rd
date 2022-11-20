\name{generateNetworkGraphPlots}
\alias{generateNetworkGraphPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Generate Multiple Plots of a Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Generate multiple plots using an igraph and corresponding node-level meta data.
}
\usage{
generateNetworkGraphPlots(
  igraph, data, print_plots = TRUE,
  plot_title = NULL, plot_subtitle = NULL,
  color_nodes_by = NULL, color_scheme = "default",
  color_legend = "auto", color_title = "auto",
  edge_width = 0.1, size_nodes_by = 0.5,
  node_size_limits = NULL, size_title = "auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{igraph}{
%%     ~~Describe \code{igraph} here~~
  An \code{igraph} object containing the edge list for the network.
}
  \item{data}{
%%     ~~Describe \code{data} here~~A data
    A data frame containing the node-level meta-data for the network, with each row corresponding to a node.
}
  \item{print_plots}{
%%     ~~Describe \code{print_plots} here~~
    A logical scalar; should visualizations be printed in the \code{R} plotting window?
}
  \item{plot_title}{
%%     ~~Describe \code{plot_title} here~~
    Passed to \code{plotNetworkGraph()}.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
    Passed to \code{plotNetworkGraph()}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Optional. Column name or number of \code{data} used to color network nodes
    in visualizations. If a vector of values is supplied, one plot will be
    generated for each value.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    Passed to \code{plotNetworkGraph}. If a vector is supplied to \code{color_nodes_by},
    this argument will optionally accept a vector of matching length; if supplied,
    corresponding values of \code{color_scheme} and \code{color_nodes_by} will
    be used together in each call to \code{plotNetworkGraph}.
}
  \item{color_legend}{
%%     ~~Describe \code{color_legend} here~~
    Passed to \code{plotNetworkGraph()}.
}
  \item{color_title}{
%%     ~~Describe \code{color_title} here~~
    Passed to \code{plotNetworkGraph}. If a vector is supplied to \code{color_nodes_by},
    this argument will optionally accept a vector of matching length; if supplied,
    corresponding values of \code{color_title} and \code{color_nodes_by} will
    be used together in each call to \code{plotNetworkGraph}.
}
  \item{edge_width}{
%%     ~~Describe \code{edge_width} here~~
    Passed to \code{plotNetworkGraph}.
}
  \item{size_nodes_by}{
%%     ~~Describe \code{size_nodes_by} here~~
    Numeric value specifying the node size in visualizations, or a
    column name or number of \code{data} for dynamically sizing the network nodes
    in visualizations.
}
  \item{node_size_limits}{
%%     ~~Describe \code{node_size_limits} here~~
    Passed to \code{plotNetworkGraph()}.
}
  \item{size_title}{
%%     ~~Describe \code{size_title} here~~
    Passed to \code{plotNetworkGraph()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  This function's primary use is to generate multiple plots of a single network, each using a different variable to color-code the nodes.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A list containing one item for each plot generated. Each plot is stored as a \code{ggraph} object.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
%% ~put references to the literature/web site here ~
    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{plotNetworkGraph}}
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
net <- generateNetworkFromSeqs(data$clone_seq,
                               drop_isolated_nodes = FALSE,
                               output_dir = NULL)

# Plot network graph
net_plot <- generateNetworkGraphPlots(
  net, color_nodes_by = c("group", "count"))
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
