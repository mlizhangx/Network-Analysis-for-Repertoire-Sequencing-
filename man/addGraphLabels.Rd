\name{labelNodes}
\alias{labelNodes}
\alias{addGraphLabels}
\title{
  Label Nodes in a Network Graph Plot
}
\description{
  Functions for annotating a graph plot to add custom labels to the nodes.
}
\usage{
labelNodes(
  net,
  node_labels,
  plots = NULL,
  size = 5,
  color = "black"
)

addGraphLabels(
  plot,
  node_labels,
  size = 5,
  color = "black"
)
}
\arguments{
  \item{net}{
  A \code{\link{list}} of network objects conforming to the output of
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
      See details.
}
  \item{plot}{
  A \code{\link[ggraph]{ggraph}} object containing the network graph plot.
}
  \item{node_labels}{
  A vector containing the node labels, where each entry is the label for a
  single node. The length should match the number of nodes in the plot.
}
  \item{plots}{
  Specifies which plots in \code{net$plots} to annotate.
  Accepts a character vector of element names or a numeric vector of element
  position indices.
  The default \code{NULL} annotates all plots.
}
  \item{size}{
  The font size of the node labels. Passed to the \code{size} argument of
  \code{\link[ggraph:geom_node_text]{geom_node_text()}}.
}
  \item{color}{
  The color of the node labels. Passed to the \code{size} argument of
  \code{\link[ggraph:geom_node_text]{geom_node_text()}}.
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

  Labels are added using
  \code{\link[ggraph:geom_node_text]{geom_node_text()}}.
}
\value{
  \code{labelNodes()} returns a copy of \code{net} with the specified plots
  annotated.

  \code{addGraphLabels()} returns
  a \code{\link[ggraph]{ggraph}} object containing the original plot annotated
  with the node labels.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
set.seed(42)
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
    seq_col = "CloneSeq"
  )

# Label each node with its receptor sequence
network <- labelNodes(network, "CloneSeq")

# View the plot
network$plots[[1]]

# Alternative method, working with the plot directly
addGraphLabels(
  network$plots[[1]],
  network$node_data$CloneSeq,
  size = 4
)
}

