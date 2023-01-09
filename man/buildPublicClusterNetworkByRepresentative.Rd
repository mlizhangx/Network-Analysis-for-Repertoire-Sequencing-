\name{buildPublicClusterNetworkByRepresentative}
\alias{buildPublicClusterNetworkByRepresentative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Build Public Cluster Network By Representative
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given cluster-level meta data for each sample's filtered public clusters,
    combine the data across samples and perform network analysis using a
    representative sequence from each cluster.
}
\usage{
buildPublicClusterNetworkByRepresentative(

  ## Input ##
  file_list =
    list.files(file.path(getwd(), "public_clusters", "cluster_meta_data")),
  input_type = "rds", data_symbols = "cdat", header = TRUE, sep = "",
  seq_col = "seq_w_max_count", count_col = "agg_count",

  ## Network ##
  dist_type = "hamming", dist_cutoff = 1,

  ## Visualization ##
  plots = TRUE, print_plots = TRUE,
  plot_title = "auto", plot_subtitle = "auto",
  color_nodes_by = c("ClusterIDPublic", "SampleID"),
  color_scheme = "turbo", color_title = c("public cluster", "sample"),
  ...,

  ## Output ##
  output_dir = file.path(getwd(), "public_clusters"),
  output_type = "rda", output_name = "PubClustByRepresentative",
  pdf_width = 12, pdf_height = 10

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} when loading the cluster-level meta data
  for each sample (one data frame per sample). Data frames should match the
  format of the cluster-level meta data returned by \code{buildRepSeqNetwork}.
}
  \item{input_type}{
%%     ~~Describe \code{input_type} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{data_symbols}{
%%     ~~Describe \code{data_symbols} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of each data frame that contains the representative sequence to be used as the basis of similarity between rows during network analysis. This column must have the same name in each sample.
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{dist_cutoff}{
%%     ~~Describe \code{dist_cutoff} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{print_plots}{
%%     ~~Describe \code{print_plots} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{plot_title}{
%%     ~~Describe \code{plot_title} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Passed to \code{generateNetworkGraphPlots}.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    Passed to \code{generateNetworkGraphPlots}.
}
  \item{color_title}{
%%     ~~Describe \code{color_title} here~~
    Passed to \code{generateNetworkGraphPlots}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
    Other arguments to \code{generateNetworkGraphPlots}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    Passed to \code{saveNetwork}.
}
  \item{output_type}{
%%     ~~Describe \code{output_type} here~~
    Passed to \code{saveNetwork}.
}
  \item{output_name}{
%%     ~~Describe \code{output_name} here~~
    Passed to \code{saveNetwork}.
}
  \item{pdf_width}{
%%     ~~Describe \code{pdf_width} here~~
    Passed to \code{saveNetwork}.
}
  \item{pdf_height}{
%%     ~~Describe \code{pdf_height} here~~
    Passed to \code{saveNetwork}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    This function is intended for performing network analysis on the public clusters
    obtained using \code{findPublicClusters()}.

    By using the cluster-level meta data as the input to \code{buildRepSeqNetwork()},
    this function treats each public cluster as a single node and performs network
    analysis with similarity based on a representative sequence for each cluster
    (e.g., the sequence with the greatest clone count). \code{buildRepSeqNetwork} is
    called without any filtering and with \code{drop_isolated_nodes = FALSE}.

    All node-level properties are automatically computed for the network, and
    cluster-level properties for the network (in which clusters of nodes represent
    clusters of public-cluster-representatives) are computed based on the properties
    of the public cluster represented by each node. See the 'value' section for a list and      description of each property.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A list of network objects as per \code{buildRepSeqNetwork()}. The data frame
    \code{cluster_data} includes the following variables:
    \item{TotalSampleLevelNodes}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the sum of the sample-network-level nodes
        contained in all such public clusters for the current row.}
    \item{TotalCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the sum of the aggregate clone count
        across all such public clusters for the current row.}
    \item{MeanOfMeanSeqLength}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the mean value of the mean sequence length
        across all such public clusters for the current row.}
    \item{MeanDegreeInPublicNet}{
        The mean network degree of all the network nodes in the cluster for
        the current row.}
    \item{MaxDegreeInPublicNet}{
        The maximum network degree of all the network nodes in the cluster for
        the current row.}
    \item{SeqWithMaxDegree}{
        The representative sequence of the network node with maximum network degree
        of all the network nodes in the cluster for the current row. If more than
        one node attains the maximum network degree, the first sequence found is
        returned.}
    \item{MaxCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this is the maximum value of the
        maximum clone count property across all such public clusters for the
        current row.}
    \item{SampleWithMaxCloneCount}{
        The sample possessing the public cluster with the maximum value
        of the maximum clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sample ID found is returned.}
    \item{SeqWithMaxCloneCount}{
        The representative sequence of the public cluster with the maximum value
        of the maximum clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sequence found is returned.}
    \item{MaxAggCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this is the maximum value of the
        aggregate clone count property across all such public clusters for the
        current row.}
    \item{SampleWithMaxAggCloneCount}{
        The sample possessing the public cluster with the maximum value
        of the aggregate clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sample ID found is returned.}
    \item{SeqWithMaxAggCloneCount}{
        The representative sequence of the public cluster with the maximum value
        of the aggregate clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sequence found is returned.}
    \item{DiameterLength}{
        See \code{?getClusterStats()}.}
    \item{Assortativity}{
        See \code{?getClusterStats()}.}
    \item{GlobalTransitivity}{
        See \code{?getClusterStats()}.}
    \item{EdgeDensity}{
        See \code{?getClusterStats()}.}
    \item{DegreeCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{ClosenessCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{EigenCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{EigenCentralityEigenvalue}{
        See \code{?getClusterStats()}.}
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{findPublicClusters}}
\code{\link{buildPublicClusterNetwork}}
}
\examples{
## Generate some data
set.seed(42)
num_samples <- 30; sample_size <- 30
groups <- sample(c(0, 1), size = num_samples, replace = TRUE)
init_chars <- c("AAAA", "AABA", "AACA", "AADA",
                "AAD", "AAA", "AABAA", "AAAAA")
init_probs_g0 <- rep(1 / length(init_chars), length(init_chars))
init_probs_g1 <- c(1, 5, 1, 1, 1, 1, 5, 1)
init_probs_g1 <- init_probs_g1 / sum(init_probs_g1)
samples <- matrix(nrow = sample_size, ncol = num_samples)
for (j in 1:num_samples) {
  init_probs <- init_probs_g0
  if (groups[[j]] == 1) { init_probs <- init_probs_g1 }
  samples[ , j] <- sample(init_chars, sample_size,
                          replace = TRUE, prob = init_probs)
}
num_edits <- 4
edit_ops <- c("insertion", "deletion", "transmutation")
edit_probs <- c(10/20, 2/20, 8/20)
chars <- c("A", "B", "C", "D")
char_probs_g0 <- c(1/2, 1/6, 1/6, 1/6)
char_probs_g1 <- c(1/3, 1/6, 1/6, 1/3)
for (k in 1:num_edits) {
  for (j in 1:num_samples) {
    char_probs <- char_probs_g0
    if ((groups[[j]] == 1)) { char_probs <- char_probs_g1 }
    for (i in 1:sample_size) {
      pos <- sample(1:nchar(samples[i, j]), size = 1)
      op <- sample(edit_ops, size = 1, prob = edit_probs)
      if (op == "insertion") {
        char <- sample(chars, size = 1, prob = char_probs)
        prefix <- substr(samples[i, j], 1, pos)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, char, suffix)
      } else if (op == "deletion") {
        prefix <- substr(samples[i, j], 1, pos - 1)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, suffix)
      } else if (op == "transmutation") {
        char <- sample(chars, size = 1, prob = char_probs)
        prefix <- substr(samples[i, j], 1, pos - 1)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, char, suffix)
      }
    }
  }
}
sample_symbols <- paste0("sample", 1:num_samples)
data_dir <- tempdir()
for (j in 1:num_samples) {
  dat <- data.frame(
    "seq" = samples[ , j], "sample_id" = j, "group" = groups[[j]],
    "count" = rbinom(sample_size, size = 1e5, prob = 0.5))
  assign(x = sample_symbols[[j]], value = dat)
  saveRDS(dat, file = file.path(data_dir, paste0(j, ".rds")))
  }

## Find public clusters
findPublicClusters(
  file_list = file.path(data_dir, paste0(1:num_samples, ".rds")),
  input_type = "rds", seq_col = "seq", count_col = "count",
  output_dir = file.path(data_dir, "filtered_samples"))

## Build public cluster network by representative sequence
input_dir <- file.path(data_dir, "filtered_samples", "cluster_meta_data")
buildPublicClusterNetworkByRepresentative(
  file_list = file.path(input_dir, list.files(input_dir)),
  output_dir = NULL)
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
