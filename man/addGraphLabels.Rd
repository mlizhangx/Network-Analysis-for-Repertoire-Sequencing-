\name{addGraphLabels}
\alias{addGraphLabels}
\title{
  Label Nodes in a Network Graph Plot
}
\description{
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
\arguments{
  \item{plot}{
  A \code{\link[ggraph]{ggraph}} object containing the network graph plot.
}
  \item{node_labels}{
  A vector containing the node labels, where each entry is the label for a single node. The length should match the number of nodes in the plot.
}
  \item{size}{
  The font size of the node labels. Passed to the \code{size} argument of \code{\link[ggraph]{geom_node_text}}.
}
  \item{color}{
  The color of the node labels. Passed to the \code{size} argument of \code{\link[ggraph]{geom_node_text}}.
}
}
\details{
  Labels are added using \code{\link[ggraph]{geom_node_text}}.
}
\value{
  A \code{\link[ggraph]{ggraph}} object containing the original plot annotated with the node labels.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
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