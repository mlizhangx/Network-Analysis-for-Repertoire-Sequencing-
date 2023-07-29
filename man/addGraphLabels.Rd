\name{addGraphLabels}
\alias{addGraphLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Label Nodes in a Network Graph Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Annotate a graph plot to add custom labels to the nodes.
}
\usage{
addGraphLabels(
  plot,
  node_labels,
  size = 5,
  color = "black"
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plot}{
%%     ~~Describe \code{plot} here~~
  A \code{\link[ggraph]{ggraph}} object containing the network graph plot.
}
  \item{node_labels}{
%%     ~~Describe \code{node_labels} here~~
  A vector containing the node labels, where each entry is the label for a single node. The length should match the number of nodes in the plot.
}
  \item{size}{
%%     ~~Describe \code{size} here~~
  The font size of the node labels. Passed to the \code{size} argument of \code{\link[ggraph]{geom_node_text}}.
}
  \item{color}{
%%     ~~Describe \code{color} here~~
  The color of the node labels. Passed to the \code{size} argument of \code{\link[ggraph]{geom_node_text}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  Labels are added using \code{\link[ggraph]{geom_node_text}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{\link[ggraph]{ggraph}} object containing the original plot annotated with the node labels.
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
toy_data <-
  simulateToyData(
    samples = 1,
    sample_size = 10,
    prefix_length = 1
  )

# Generate network
network <-
  buildRepSeqNetwork(
    toy_data,
    seq_col = "CloneSeq",
    print_plots = FALSE,
    output_dir = NULL
  )

# Label each node with its receptor sequence
network$plots[[1]] <-
  addGraphLabels(
    network$plots[[1]],
    network$node_data$CloneSeq,
    size = 4
  )

# View the plot
network$plots[[1]]
}