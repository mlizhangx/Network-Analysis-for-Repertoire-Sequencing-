\name{buildRepSeqNetwork}
\alias{buildRepSeqNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Immune Repertoire Network By Sequence Similarity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Builds the network graph for an immune repertoire based on sequence
    similarity, computes specified network properties and generates customized
    visualizations.
}
\usage{
buildRepSeqNetwork(

  ## Input ##
  data, seq_col, count_col = NULL,
  subset_cols = NULL,
  min_seq_length = 3, drop_matches = NULL,

  ## Network ##
  dist_type = "hamming", dist_cutoff = 1,
  drop_isolated_nodes = TRUE,
  node_stats = FALSE, stats_to_include = node_stat_settings(),
  cluster_stats = FALSE,

  ## Visualization ##
  plots = TRUE, print_plots = TRUE,
  plot_title = "auto", plot_subtitle = "auto",
  color_nodes_by = "auto", ...,

  ## Output ##
  output_dir = getwd(), output_type = "individual",
  output_name = "MyRepSeqNetwork",
  pdf_width = 12, pdf_height = 10

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
    A data frame containing the immune repertoire sequencing data, with variables indexed by column and observations (e.g., clones or cells) indexed by row.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of \code{data} containing the receptor sequences to be used as the basis of similarity between rows. Also accepts a vector of length 2 specifying distinct sequence columns (e.g., alpha, beta chains); then two rows are similar only if both types of sequences are similar.
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
    Optional column name or number of \code{data} containing a measure of abundance, e.g., clone count. Passed to \code{getClusterStats}; only relevant if \code{cluster_stats = TRUE}.
}
  \item{subset_cols}{
%%     ~~Describe \code{other_cols} here~~
    Optional vector of column names or numbers of \code{data}; if supplied, only these plus other relevant columns are included in the output. If \code{NULL} (the default), all columns of \code{data} are included.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    A numeric scalar, or \code{NULL}. Observations whose receptor sequences have fewer than \code{min_seq_length} characters are removed prior to network analysis.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_chars} here~~
    An optional regular expression or character string. Observations whose receptor sequences return a match are removed prior to network analysis.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    The type of function to use as a measure of similarity between two receptor
    sequences. Valid options are \code{"hamming"} (the default), \code{"levenshtein"}
    and \code{"euclidean_on_atchley"} (only applicable to TCR CDR3 amino acid sequences).
}
  \item{dist_cutoff}{
%%     ~~Describe \code{dist_cutoff} here~~
    A nonnegative scalar specifying the maximum distance threshold for
    similarity between receptor sequences.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
    A logical scalar; should observations whose receptor sequences are not
    similar to any other sequences be dropped from the network?
}
  \item{node_stats}{
%%     ~~Describe \code{node_stats} here~~
    A logical scalar; should node-level network properties be computed?
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
    A list returned by \code{node_stat_settings()}, specifying the node-level
    properties to include. Also accepts the values \code{"all"} and
    \code{"cluster_id_only"}. Only relevant if \code{node_stats = TRUE}.
}
  \item{cluster_stats}{
%%     ~~Describe \code{cluster_stats} here~~
    A logical scalar; should cluster-level network properties be computed?
}
  \item{plots}{
    A logical scalar; should plot(s) of the network graph be generated?
}
  \item{print_plots}{
%%     ~~Describe \code{print_plots} here~~
    A logical scalar; should visualizations be printed in the \code{R} plotting window?
}
  \item{plot_title}{
%%     ~~Describe \code{plot_title} here~~
    A character string (or \code{NULL} value) to be used as the title in
    visualizations. The default value \code{"auto"} generates the title
    dynamically.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
    A character string (or \code{NULL} value) to be used as the subtitle in
    visualizations. The default value \code{"auto"} generates the title
    dynamically.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Optional. Column name or number of \code{data} used to color network nodes
    in visualizations. If a vector of values is supplied, one plot will be
    generated for each value. The value \code{"auto"} attempts to use one of several        variables (e.g., cluster ID) to color the nodes, depending on what is available.
}
  \item{...}{
  Other arguments to \code{generateNetworkGraphPlots}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
  The output directory; if \code{NULL}, output will not be written to file.
}
  \item{output_type}{
%%     ~~Describe \code{output_dir} here~~
  A character string specifying the file format to use when writing output to file. Default \code{"individual"} saves each item as a separate, uncompressed file, with data frames saved in csv format. \code{"rda"} and \code{"rds"} save the output list as a rda and rds file, respectively. For all output types, any plots generated will also be saved in a pdf file.
}
  \item{output_name}{
%%     ~~Describe \code{output_dir} here~~
  A character string to be used as a common filename prefix for any files saved.
}
  \item{pdf_width}{
%%     ~~Describe \code{plot_width} here~~
    Passed to the \code{width} argument of \code{grDevices::pdf()}.
}
  \item{pdf_height}{
%%     ~~Describe \code{plot_height} here~~
    Passed to the \code{height} argument of \code{grDevices::pdf()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    To build the immune repertoire network, each TCR/BCR clone (bulk data) or
    cell (single-cell data) is modeled as a network node (corresponding to a single
    row of the input data). Two nodes are considered adjacent (share an edge) if
    their receptor sequences are sufficiently similar.

    Both nucleotide and amino acid sequences are supported. Sequence similarity is
    based on either the Hamming distance or Levenshtein (edit) distance. For TCR CDR3
    amino acid sequences, an alternative measure of similarity is also available via
    the argument \code{dist_type = "euclidean_on_atchley"}; this represents sequences as
    30-dimensional numeric vectors according to the Atchley factors of their amino acids
    (which encode biological properties) using a trained encoder. Similarity is then based on the
    Euclidean distance between these representations. Selecting a distance
    threshold is challenging in this setting, however, as particular values lack
    the clear interpretations enjoyed by the Hamming and Levenshtein distances.

    The graph adjacency matrix for the network is computed and used to generate
    a \code{igraph} network object. Network properties are computed largely
    through the use of \code{igraph} functions. The visualizations are generated with
    \code{ggraph}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A list containing some or all of the following items:

    \item{igraph}{An \code{igraph} object containing the edge list for the network.}
    \item{adjacency_matrix}{The network graph adjacency matrix, stored as a sparse matrix of class \code{dgCMatrix} from the \code{Matrix} package.}
    \item{node_data}{ A data
    frame containing the node-level meta-data for the network. This data frame contains
    all column names of \code{data} unless \code{subset_cols} is non-null, in which
    case only the column names of \code{data} specified by \code{subset_cols}, along
    with other relevant ones (e.g., \code{seq_col}), will be present in the returned
    data frame. The data frame will additionally contain columns for any node-level
    network properties computed (as specified through the \code{node_stats} and
    \code{stats_to_include} arguments). The data frame will contain one row for
    each node that remains in the network after all filtering has occurred,
    including any filtering based on \code{drop_isolated_nodes}. The row names
    of the original input data will be preserved. }
    \item{cluster_data}{A data frame containing the cluster-level properties, with one row per cluster. Only included if \code{cluster_stats = TRUE}.}
    \item{plots}{A list containing the graph plot(s) as \code{ggraph} objects. Only included if \code{plots = TRUE}.}

    The output of \code{buildRepSeqNetwork} is returned invisibly (it will not be printed to console when the output is not assigned).
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
# Generate some data
set.seed(42)
sample_size <- 200
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
data <- data.frame(seqs = clone_seqs)
net <- buildRepSeqNetwork(data, "seqs", output_dir = NULL)

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
