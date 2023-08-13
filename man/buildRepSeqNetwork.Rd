\name{buildRepSeqNetwork}
\alias{buildRepSeqNetwork}
\title{
    Network Analysis of Immune Repertoire
}
\description{
    Given Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data, builds
    the network graph for the immune repertoire based on sequence similarity,
    computes specified network properties and generates customized visualizations.
}
\usage{
buildRepSeqNetwork(

  ## Input ##
  data,
  seq_col,
  count_col = NULL,
  subset_cols = NULL,
  min_seq_length = 3,
  drop_matches = NULL,

  ## Network ##
  dist_type = "hamming",
  dist_cutoff = 1,
  drop_isolated_nodes = TRUE,
  node_stats = FALSE,
  stats_to_include = chooseNodeStats(),
  cluster_stats = FALSE,
  cluster_fun = cluster_fast_greedy,

  ## Visualization ##
  plots = TRUE,
  print_plots = TRUE,
  plot_title = "auto",
  plot_subtitle = "auto",
  color_nodes_by = "auto",
  ...,

  ## Output ##
  output_dir = getwd(),
  output_type = "individual",
  output_name = "MyRepSeqNetwork",
  pdf_width = 12,
  pdf_height = 10

)
}
\arguments{
  \item{data}{
    A data frame containing the AIRR-Seq data, with variables indexed by column
    and observations (e.g., clones or cells) indexed by row.
}
  \item{seq_col}{
    Specifies the column of \code{data} containing the receptor sequences to be
    used as the basis of similarity between rows. Accepts either the column name
    or column index. Alternatively, a vector of length 2 specifying two columns
    with distinct sequences (e.g., alpha chain and beta chain), in which case
    similarity between rows depends on both types of sequences being similar
    (e.g., similar alpha chain sequences and similar beta chain sequences).
}
  \item{count_col}{
    Optional. Specifies the column of \code{data} containing a measure of abundance,
    e.g., clone count or unique molecular identifier (UMI) count. Accepts either
    the column name or column index. Passed to \code{\link{getClusterStats}}; only
    relevant if \code{cluster_stats = TRUE}.
}
  \item{subset_cols}{
    Specifies which columns of the AIRR-Seq data are included in the output.
    Accepts a vector of column names or a vector of column indices. The default
    \code{NULL} includes all columns. The receptor sequence column is always
    included regardless of this argument's value.
    Passed to \code{\link{filterInputData}}.
}
  \item{min_seq_length}{
    A numeric scalar, or \code{NULL}. Observations whose receptor sequences have
    fewer than \code{min_seq_length} characters are removed prior to network analysis.
}
  \item{drop_matches}{
    Optional. Passed to \code{\link{filterInputData}}. Accepts a regular expression
    (see \code{\link[base]{regex}}). Checks receptor sequences for a pattern match
    using \code{\link{grep}}. Those returning a match are removed prior to network
    analysis.
}
  \item{dist_type}{
    Specifies the function used to quantify the similarity between sequences.
    The similarity between two sequences determines the pairwise distance between
    their respective nodes in the network graph, with greater similarity corresponding
    to shorter distance. Valid options are \code{"hamming"} (the default), which
    uses \code{\link{hamDistBounded}}, and \code{"levenshtein"}, which uses
    \code{\link{levDistBounded}}.
}
  \item{dist_cutoff}{
    A nonnegative scalar. Specifies the maximum pairwise distance (based on
    \code{dist_type}) for an edge connection to exist between two nodes. Pairs of
    nodes whose distance is less than or equal to this value will be joined by an
    edge connection in the network graph. Controls the stringency of the network
    construction and affects the number and density of edges in the network. A lower
    cutoff value requires greater similarity between sequences in order for their
    respective nodes to be joined by an edge connection. A value of \code{0}
    requires two sequences to be identical in order for their nodes to be joined
    by an edge.
}
  \item{drop_isolated_nodes}{
    A logical scalar. When \code{TRUE}, removes each node that is not joined by an
    edge connection to any other node in the network graph.
}
  \item{node_stats}{
    A logical scalar. Specifies whether node-level network properties are computed.
}
  \item{stats_to_include}{
    A named logical vector returned by \code{\link{chooseNodeStats}} or
    \code{\link{exclusiveNodeStats}}. Specifies the node-level network properties
    to compute. Also accepts the values \code{"all"} and \code{"cluster_id_only"}.
    Only relevant if \code{node_stats = TRUE}.
}
  \item{cluster_stats}{
    A logical scalar. Specifies whether to compute cluster-level network properties.
}
  \item{cluster_fun}{
    Passed to \code{\link{addClusterMembership}}. Specifies the clustering algorithm
    used when cluster analysis is performed. Cluster analysis is performed when
    \code{cluster_stats = TRUE} or when \code{node_stats = TRUE} with the
    \code{cluster_id} property enabled via the \code{stats_to_include} argument.
}
  \item{plots}{
    A logical scalar. Specifies whether to generate plots of the network graph.
}
  \item{print_plots}{
    A logical scalar. If \code{plots = TRUE}, specifies whether the plots should
    be printed to the R plotting window.
}
  \item{plot_title}{
    A character string or \code{NULL}. If \code{plots = TRUE}, this is the title
    used for each plot. The default value \code{"auto"} generates the title based
    on the value of the \code{output_name} argument.
}
  \item{plot_subtitle}{
    A character string or \code{NULL}. If \code{plots = TRUE}, this is the subtitle
    used for each plot. The default value \code{"auto"} generates a subtitle based
    on the values of the \code{dist_type} and \code{dist_cutoff} arguments.
}
  \item{color_nodes_by}{
    Optional. Specifies a variable to be used as metadata for coloring the nodes
    in the network graph plot. Accepts a character string. This can be a column
    name of \code{data} or (if \code{node_stats = TRUE}) the name of a computed
    node-level network property (based on \code{stats_to_include}). Also accepts
    a character vector specifying multiple variables, in which case one plot will
    be generated for each variable. The default value \code{"auto"} attempts to use
    one of several potential variables to color the nodes, depending on what is
    available. A value of \code{NULL} leaves the nodes uncolored.
}
  \item{...}{
  Other arguments to \code{\link{generateNetworkGraphPlots}}.
}
  \item{output_dir}{
  A file path specifying the directory for saving the output. The directory will
  be created if it does not exist. If \code{NULL}, output will be returned but
  not saved.
}
  \item{output_type}{
  A character string specifying the file format to use when saving the output.
  The default value \code{"individual"} saves each element of the returned list
  as an individual uncompressed file, with data frames saved in csv format. For
  better compression, the values \code{"rda"} and \code{"rds"} save the returned
  list as a single file using the rda and rds format, respectively (in the former
  case, the list will be named \code{net} within the rda file). Regardless of the
  argument value, any plots generated will saved to a pdf file containing one plot
  per page.
}
  \item{output_name}{
  A character string. All files saved will have file names beginning with this value.
}
  \item{pdf_width}{
    Sets the width of each plot when writing to pdf.
    Passed to \code{\link{saveNetwork}}.
}
  \item{pdf_height}{
    Sets the height of each plot when writing to pdf.
    Passed to \code{\link{saveNetwork}}.
}
}
\details{
To construct the immune repertoire network, each TCR/BCR clone (bulk data)
or cell (single-cell data) is modeled as a node in the network graph,
corresponding to a single row of the AIRR-Seq data. For each node, the
corresponding receptor sequence is considered. Both nucleotide and amino acid
sequences are supported for this purpose. The receptor sequence is used as the
basis of similarity and distance between nodes in the network.

Similarity between sequences is measured using either the Hamming distance
or Levenshtein (edit) distance. The similarity determines the pairwise distance
between nodes in the network graph. The more similar two sequences are, the
shorter the distance between their respective nodes. Two nodes in the graph
are joined by an edge if the distance between them is sufficiently small,
i.e., if their receptor sequences are sufficiently similar.

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html}{buildRepSeqNetwork package vignette}
for more details. The vignette can be accessed offline using
\code{vignette("buildRepSeqNetwork")}.
}
\value{
    If the constructed network contains no nodes, the function will return
    \code{NULL}, invisibly, with a warning. Otherwise, the function invisibly
    returns a list containing the following items:

    \item{igraph}{An object of class \code{\link[igraph]{igraph}} containing the
    list of nodes and edges for the network graph.
    }
    \item{adjacency_matrix}{The network graph adjacency matrix, stored as a sparse
    matrix of class \code{dgCMatrix} from the \code{Matrix} package.
    See \code{\link[Matrix]{dgCMatrix-class}}.
    }
    \item{node_data}{ A data frame containing containing metadata for the network
    nodes, where each row corresponds to a node in the network graph. This data
    frame contains all variables from \code{data} (unless otherwise specified via
    \code{subset_cols}) in addition to the computed node-level network properties
    if \code{node_stats = TRUE}. Each row's name is the name of the corresponding
    row from \code{data}.
    }
    \item{cluster_data}{A data frame containing network properties for the clusters,
    where each row corresponds to a cluster in the network graph. Only included if
    \code{cluster_stats = TRUE}.
    }
    \item{plots}{A list containing each plot generated for the network graph.
    Each plot is an object of class \code{\link[ggraph]{ggraph}}. Only included
    if \code{plots = TRUE}.
    }
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html}{buildRepSeqNetwork vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{

toy_data <- simulateToyData()

network <- buildRepSeqNetwork(
  toy_data, seq_col = "CloneSeq",
  node_stats = TRUE, cluster_stats = TRUE,
  color_nodes_by = "transitivity", color_scheme = "plasma-1",
  size_nodes_by = "degree", node_size_limits = c(0.5, 1.5),
  output_dir = NULL)


}