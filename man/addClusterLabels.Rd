\name{labelClusters}
\alias{labelClusters}
\alias{addClusterLabels}
\title{
  Label Clusters in a Network Graph Plot
}
\description{
  Functions for labeling the clusters in network graph plots with their cluster
  IDs.
  The user can specify a cluster-level property by which to rank the clusters,
  labeling only those clusters above a specified rank.
}
\usage{
labelClusters(
  net,
  plots = NULL,
  top_n_clusters = 20,
  cluster_id_col = "cluster_id",
  criterion = "node_count",
  size = 5, color = "black",
  greatest_values = TRUE
)

addClusterLabels(
  plot,
  net,
  top_n_clusters = 20,
  cluster_id_col = "cluster_id",
  criterion = "node_count",
  size = 5,
  color = "black",
  greatest_values = TRUE
)
}
\arguments{
  \item{net}{
  A \code{\link{list}} of network objects conforming to the output of
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
      See details.
}
  \item{plots}{
  Specifies which plots in \code{net$plots} to annotate.
  Accepts a character vector of element names or a numeric vector of element
  position indices.
  The default \code{NULL} annotates all plots.
}
  \item{plot}{
  A \code{\link[ggraph]{ggraph}} object containing the network graph plot.
}
  \item{top_n_clusters}{
  A positive integer specifying the number of clusters to label. Those with the
  highest rank according to the \code{criterion} argument will be labeled.
}

  \item{cluster_id_col}{
  Specifies the column of \code{net$node_data} containing
  the variable for cluster membership.
  Accepts a character string containing the column name.
}
  \item{criterion}{
  Can be used to specify a cluster-level network property by which to rank the
  clusters. Non-default values are ignored unless \code{net$cluster_data} exists and
  corresponds to the cluster membership variable specified by \code{cluster_id_col}.
  Accepts a character string containing a column name of \code{net$cluster_data}.
  The property must be quantitative for the ranking to be meaningful.
  By default, clusters are ranked by node count, which is computed based on the
  cluster membership values if necessary.
}
  \item{size}{
  The font size of the cluster ID labels. Passed to the \code{size} argument of
  \code{\link[ggraph:geom_node_text]{geom_node_text()}}.
}
  \item{color}{
  The color of the cluster ID labels. Passed to the \code{color} argument of
  \code{\link[ggraph:geom_node_text]{geom_node_text()}}.
}
  \item{greatest_values}{
  Logical. Controls whether clusters are ranked according to the greatest or
  least values of the property specified by the \code{criterion} argument. If
  \code{TRUE}, clusters with greater values will be ranked above those with
  lower values, thereby receiving a higher priority to be labeled.
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
}
\value{
  \code{labelClusters()} returns a copy of \code{net} with the specified plots
  annotated.

  \code{addClusterLabels()} returns an annotated copy of \code{plot}.
}
\seealso{
\code{\link[=addClusterMembership]{addClusterMembership()}},
\code{\link[=getClusterStats]{getClusterStats()}},
\code{\link[=generateNetworkGraphPlots]{generateNetworkGraphPlots()}}
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
toy_data <- simulateToyData()

network <- buildRepSeqNetwork(
  toy_data, "CloneSeq",
  cluster_stats = TRUE,
  color_nodes_by = "cluster_id",
  color_scheme = "turbo",
  color_legend = FALSE,
  plot_title = NULL,
  plot_subtitle = NULL,
  size_nodes_by = 1
)

network <- labelClusters(network)

network$plots$cluster_id
}
