\name{buildPublicClusterNetwork}
\alias{buildPublicClusterNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Build Network of Public Clusters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given node-level meta data for each sample's filtered public clusters,
    combine the data across samples and perform network analysis.
}
\usage{
buildPublicClusterNetwork(

  ## Input ##
  file_list =
    list.files(file.path(getwd(), "public_clusters", "node_meta_data")),
  input_type = "rds", data_symbols = "ndat", header = TRUE, sep = "",
  seq_col,

  ## Network ##
  drop_isolated_nodes = FALSE,
  node_stats = TRUE, stats_to_include = "all", cluster_stats = TRUE,

  ## Visualization ##
  color_nodes_by = c("ClusterIDPublic", "SampleID"),
  color_scheme = "turbo", color_title = c("public cluster", "sample"),

  ## Output ##
  output_dir = file.path(getwd(), "public_clusters"),
  output_name = "PublicClusterNetwork",

  ...

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} when loading the node-level meta data
  for each sample (one data frame per sample).
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
    The column name or number of each sample's data frame that contains the receptor sequences to be used as the basis of similarity between rows during network analysis. This column must have the same name in each sample.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{node_stats}{
%%     ~~Describe \code{node_stats} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{cluster_stats}{
%%     ~~Describe \code{cluster_stats} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_title}{
%%     ~~Describe \code{color_title} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{output_name}{
%%     ~~Describe \code{output_name} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
    Other arguments to \code{buildRepSeqNetwork}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    The input data is intended to be obtained using \code{findPublicClusters()}
    on the unfiltered RepSeq data for each sample.

    Any node-level properties computed are renamed to reflect their association
    to the public network. The names used are
    \code{ClusterIDPublic}, \code{PublicNetworkDegree},
    \code{PublicTransitivity}, \code{PublicCloseness},
    \code{PublicCentralityByCloseness}, \code{PublicEigenCentrality},
    \code{PublicCentralityByEigen}, \code{PublicBetweenness},
    \code{PublicCentralityByBetweenness}, \code{PublicAuthorityScore},
    \code{PublicCoreness}, \code{PublicPageRank}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A list of network objects as returned by \code{buildRepSeqNetwork()}.
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
\code{\link{buildPublicClusterNetworkByRepresentative}}
}
\examples{
## Generate some data
set.seed(42)
num_samples <- 30; sample_size <- 30
groups <- sample(c(0, 1), size = num_samples, replace = TRUE)
init_chars <- c("AAAA", "AABA", "AACA", "AADA", "AAD", "AAA", "AABAA", "AAAAA")
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

## Build public cluster network
input_dir <- file.path(data_dir, "filtered_samples", "node_meta_data")
buildPublicClusterNetwork(
  file_list = file.path(input_dir, list.files(input_dir)),
  seq_col = "seq", count_col = "count", output_dir = NULL)
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
