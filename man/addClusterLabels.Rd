\name{addClusterLabels}
\alias{addClusterLabels}
\title{
  Label Clusters in a Network Graph Plot
}
\description{
  Given a graph plot and metadata for a network where clustering analysis has been performed, labels clusters in the plot with their cluster IDs. The user can specify a cluster-level property by which to rank the clusters, labeling only those clusters above a specified rank.
}
\usage{
addClusterLabels(
  plot,
  net,
  top_n_clusters = 20,
  cluster_id_col = "cluster_id",
  criterion = "node_count",
  size = 5,
  color = "black",
  greatest_values = TRUE)
}
\arguments{
  \item{plot}{
  A \code{\link[ggraph]{ggraph}} object containing the network graph plot.
}
  \item{net}{
  A list containing data frames named \code{node_data} and \code{cluster_data} with node-level and cluster-level network properties, such as that returned by \code{\link{buildRepSeqNetwork}} when called with \code{cluster_stats = TRUE}.
}
  \item{top_n_clusters}{
  A positive integer specifying the number of clusters to label. Those with the highest rank according to the \code{criterion} argument will be labeled.
}
  \item{cluster_id_col}{
  The column name of the cluster ID variable within \code{net$node_data}. Note that for \code{net$cluster_data}, the variable is always assumed to be \code{cluster_id}.
}
  \item{criterion}{
  A column name of \code{net$cluster_data} specifying a cluster-level network property by which to rank the clusters. The property should be quantitative for the ranking to be meaningful. See \code{\link{getClusterStats}} for a list of properties. The default is \code{"node_size"}.
}
  \item{size}{
  The font size of the cluster ID labels. Passed to the \code{size} argument of \code{\link[ggraph]{geom_node_text}}.
}
  \item{color}{
  The color of the cluster ID labels. Passed to the \code{color} argument of \code{\link[ggraph]{geom_node_text}}.
}
  \item{greatest_values}{
  Logical. Controls whether clusters are ranked according to the greatest or least values of the property specified by the \code{criterion} argument. If \code{TRUE}, clusters with greater values will be ranked above those with lower values, thereby receiving a higher priority to be labeled.
}
}
\value{
  A \code{\link[ggraph]{ggraph}} object containing the original plot annotated with the cluster ID labels.
}
\seealso{
\code{\link{addClusterMembership}},
\code{\link{getClusterStats}},
\code{\link{generateNetworkGraphPlots}}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData()

network <- buildRepSeqNetwork(
  toy_data,
  seq_col = "CloneSeq",
  node_stats = TRUE,
  cluster_stats = TRUE,
  color_nodes_by = "transitivity",
  color_scheme = "plasma-1",
  size_nodes_by = "degree",
  node_size_limits = c(0.5, 1.5),
  print_plots = FALSE,
  output_dir = NULL
)

# Label two largest clusters
addClusterLabels(
  plot = network$plots$transitivity,
  net = network,
  top_n_clusters = 2,
  criterion = "node_count" # (the default)
)
}
