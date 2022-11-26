\name{kmeansAtchley}
\alias{kmeansAtchley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Analyze TCR Samples Using Numerical Encoder and K-Means Clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given samples of TCR RepSeq data, convert the CDR3 amino acid
    sequences into numeric vectors using deep learning, perform \eqn{K}-means
    clustering, profile how each sample's unique TCR sequences are distributed
    among the clusters, and compare these profiles across samples.
}
\usage{
kmeansAtchley(
  data,
  amino_col = "AminoAcidSeq",
  sample_col = "SampleID",
  group_col, k = 100,
  pdf_width = 15,
  pdf_height = 15,
  margin_cluster_heatmap = 25,
  margin_corr_heatmap = 15,
  use_viridis = FALSE,
  output_dir = getwd(),
  file_cluster_heatmap =
    "atchley_kmeans_cluster_relative_size_profiles_by_sample.pdf",
  file_corr_heatmap =
    "atchley_kmeans_corr_in_cluster_size_profile_between_samples.pdf",
  return_output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
    A data frame containing the combined immune repertoire sequencing data for all samples,
    with variables indexed by column and observations indexed by row.
}
  \item{amino_col}{
%%     ~~Describe \code{amino_col} here~~
    The column of \code{data} containing the CDR3 amino acid sequences.
}
\item{sample_col}{
    The column of \code{data} containing the sample IDs.
}
\item{group_col}{
    The column of \code{data} containing the group IDs.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
    The number of clusters for \eqn{K}-means clustering.
}
  \item{pdf_width}{
%%     ~~Describe \code{plot_width} here~~
    Passed to the \code{width} argument of \code{grDevices::pdf()}.
}
  \item{pdf_height}{
%%     ~~Describe \code{plot_height} here~~
    Passed to the \code{height} argument of \code{grDevices::pdf()}.
}
  \item{margin_cluster_heatmap}{
%%     ~~Describe \code{margin_cluster_heatmap} here~~
    The margin width of column names and row names in the heatmap for each
    sample's cluster profile.
}
  \item{margin_corr_heatmap}{
%%     ~~Describe \code{margin_corr_heatmap} here~~
    The margin width of column names and row names in the heatmap for the
    correlation in cluster profiles between samples.
}
  \item{use_viridis}{
%%     ~~Describe \code{use_viridis} here~~
    A logical scalar indicating whether to use the \code{cividis}
    color scale from the \code{viridis} package instead of the \code{RdBu}
    color scale from \code{RColorBrewer}. This yields more readable results when
    printing in black-and-white, and is more robust against color
    vision deficiency.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    The directory in which to save pdf files of the two heatmaps produced.
}
  \item{file_cluster_heatmap}{
%%     ~~Describe \code{file_cluster_heatmap} here~~
    The pdf filename for the heatmap for each
    sample's cluster profile.
}
  \item{file_corr_heatmap}{
%%     ~~Describe \code{file_corr_heatmap} here~~
    The pdf filename for the heatmap for correlation in cluster profiles between samples.
}
  \item{return_output}{
%%     ~~Describe \code{return_output} here~~
    A logical scalar indicating whether to return the numeric vectors encoding
    the TCR sequences, the K-means cluster membership table, and each sample's
    profile of TCR representation across clusters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    Each unique TCR sequence is encoded as a 30-dimensional numeric vector
    using \code{encodeTCRSeqsByAtchleyFactor()}, and \eqn{K}-means clustering is
    performed on the encoded TCR sequences from all samples, resulting in
    \code{k} clusters.

    For each sample, its TCR representation profile across the \code{k} clusters
    is computed as a length-\code{k} vector whose \eqn{i}th element is the
    fraction of the sample's unique TCR sequences that belong to the \eqn{i}th
    cluster.

    Two heatmaps are generated and saved to file: the first shows the TCR
    representation profile values across samples and clusters; the second shows
    the correlation in TCR represenation profiles between pairs of samples.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If \code{return_output = TRUE}, then a list containing the following items:
    \item{kmeans_cluster_ids}{
    A data frame with two variables, \code{cdr3} and \code{kmeanClusterID},
    containing the unique TCR sequences and the cluster to which each belongs.
    }
    \item{encoded_values}{
    A matrix returned by \code{encodeTCRSeqsByAtchleyFactor} containing the
    numerically encoded values of the TCR sequences.
    }
    \item{cluster_TCR_profiles}{
    A matrix with \code{k} rows and one column per sample; each column is the
    TCR representation profile for a sample, with each value recording the
    fraction of that sample's unique TCR sequences belonging to the cluster for
    the corresponding row (thus each column sums to 1).
    }
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
\dontrun{
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

files <- file.path(data_dir, paste0(1:num_samples, ".rds"))

## Find receptor sequences associated with comparison group
associated_seqs <- findAssociatedSeqs(
  file_list = files, input_type = "rds",
  group_ids = groups,
  min_seq_length = NULL,
  drop_matches = NULL,
  groups = c(0, 1),
  seq_col = "seq",
  outfile = NULL)

## Gather clones in public neighborhoods of each associated sequence
findAssociatedClones(
  file_list = files, input_type = "rds",
  group_ids = groups,
  seq_col = "seq",
  dist_type = "levenshtein",
  assoc_seqs = associated_seqs$ReceptorSeq,
  output_dir = data_dir)

## Combine neighborhoods and perform clustering/network analysis
all_clusters <- buildAssociatedClusterNetwork(
  file_list = file.path(data_dir,
                        paste0(associated_seqs$ReceptorSeq, ".csv")),
  seq_col = "seq", dist_type = "levenshtein", size_nodes_by = 1.5,
  output_dir = NULL)

## Analysis by numeric encoding and K-means clustering
store_results <- kmeansAtchley(
  data = all_clusters$node_data,
  amino_col = "seq", sample_col = "SampleID", group_col = "GroupID",
  k = 5, output_dir = data_dir, return_output = TRUE)
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
