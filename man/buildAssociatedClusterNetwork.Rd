\name{buildAssociatedClusterNetwork}
\alias{buildAssociatedClusterNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Combine Associated Sequence Neighborhoods and Build Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    After running \code{findAssociatedClones}, this function is used to combine
    some or all of the resulting neighborhoods into a single network in which
    clustering and network analysis are performed.
}
\usage{
buildAssociatedClusterNetwork(
    file_list, input_type = "csv",
    data_symbols = NULL,
    header = TRUE, sep = ",",
    seq_col,
    min_seq_length = NULL, drop_matches = NULL,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE,
    stats_to_include =
      node_stat_settings(cluster_id = TRUE),
    cluster_stats = TRUE,
    color_nodes_by = "cluster_id",
    output_name = "AssociatedClusterNetwork",
    ...
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList}.
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
  Passed to \code{buildRepSeqNetwork}.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
  Passed to \code{buildRepSeqNetwork}.
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
    Essentially a wrapper to \code{loadDataFromFileList} and
    \code{buildRepSeqNetwork}, with default argument values tailored to the
    associated cluster workflow.
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
\code{\link{findAssociatedSeqs}}
\code{\link{findAssociatedClones}}
}
\examples{
# Generate some data
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
  min_seq_length = NULL,
  dist_type = "levenshtein",
  assoc_seqs = associated_seqs$ReceptorSeq,
  output_dir = data_dir)

## Combine neighborhoods and perform clustering/network analysis
all_clusters <- buildAssociatedClusterNetwork(
  file_list = file.path(data_dir,
                        paste0(associated_seqs$ReceptorSeq, ".csv")),
  seq_col = "seq", dist_type = "levenshtein", size_nodes_by = 1.5,
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
