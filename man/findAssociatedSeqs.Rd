\name{findAssociatedSeqs}
\alias{findAssociatedSeqs}
\alias{findAssociatedSeqs2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Find Receptor Sequences Associated to a Comparison Group
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple samples of RepSeq data from two groups, search for sequences
    associated to the comparison group based on Fisher's exact test P-value.

    \code{findAssociatedSeqs} is for use with separate files and data frames per
    sample, while \code{findAssociatedSeqs2} is for use with a single data frame
    containing all samples.
}
\usage{
findAssociatedSeqs(
  ## Input ##
  file_list, input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids = 1:length(file_list),
  subject_ids = sample_ids,
  group_ids, groups = c("group0", "group1"),
  seq_col, freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7, drop_matches = "[*|_]",
  min_sample_membership = 5, pval_cutoff = 0.05,

  ## Output ##
  outfile = "associated_seqs.csv")


findAssociatedSeqs2(
  ## Input ##
  data, seq_col, sample_col,
  subject_col = sample_col,
  group_col, groups = c("group0", "group1"),
  freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7, drop_matches = "[*|_]",
  min_sample_membership = 5, pval_cutoff = 0.05,

  ## Ouptut ##
  outfile = "associated_seqs.csv")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} (one data frame per sample).
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
  \item{sample_ids}{
%%     ~~Describe \code{sample_ids} here~~
  A character or numeric vector of sample IDs, whose length matches that of \code{file_list}.
}
  \item{subject_ids}{
%%     ~~Describe \code{subject_ids} here~~
  A character or numeric vector of subject IDs, whose length matches that of \code{file_list}.
  This can be disregarded if each sample comes from a distinct subject.
}
  \item{group_ids}{
%%     ~~Describe \code{group_ids} here~~
  A character or numeric vector of group IDs, which contains exactly two unique
  values and whose length matches that of \code{file_list}.
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
  A vector of length two having the same type as \code{group_ids} and containing the
  two unique values; the first and second elements identify the reference
  and comparison groups, respectively.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of each sample's data frame that contains the receptor sequences to associate with the comparison group.
}
  \item{freq_col}{
%%     ~~Describe \code{freq_col} here~~
    Optional. The column name or number of each sample's data frame that contains the clone frequency/clone fraction (i.e., measure of clonal abundance normalized based on the total count in the sample). This column must have the same name in each sample. If supplied, the maximum clone frequency of each associated sequence across all samples will be recorded in the returned output.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    Passed to \code{filterInputData()}.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{filterInputData()}.
}
  \item{min_sample_membership}{
%%     ~~Describe \code{min_sample_membership} here~~
    Only sequences that appear in at least this many samples will be considered.
}
  \item{pval_cutoff}{
%%     ~~Describe \code{pval_cutoff} here~~
    Only sequences with Fisher's exact test $P$-value below this cutoff will be
    returned.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
    A file path used to save the data frame containing the associated sequences
    (using \code{write.csv}).
}
\item{data}{
    A data frame containing the combined immune repertoire sequencing data for all samples,
    with variables indexed by column and observations indexed by row.
}
\item{sample_col}{
    The column of \code{data} containing the sample IDs.
}
\item{subject_col}{
    The column of \code{data} containing the subject IDs.
}
\item{group_col}{
    The column of \code{data} containing the group IDs.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The data is first filtered by minimum sequence length and sequence content
    if applicable, then is filtered based on minimum sample membership, keeping
    only those receptor sequences that appear in at least \code{min_sample_membership}
    samples.

    For each remaining sequence, Fisher's exact test is performed using
    \code{stats::fisher.test()}, based on the number of subjects possessing the
    sequence in each of the reference and comparison groups. Only sequences
    with a \eqn{P}-value below \code{pval_cutoff} will be retained. The resulting
    sequences are sorted by \eqn{P}-value in increasing order and returned
    along with some basic meta-data.

    The returned ouput is intended for downstream use with
    \code{findAssociatedClones} and \code{buildAssociatedClusterNetwork}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the following variables:
    \item{ReceptorSeq}{The receptor sequences associated with the comparison
    group, based on the specified criteria.}
    \item{shared_by_n_samples}{The number of samples in which each sequence
    appears.}
    \item{fisher_pvalue}{The P-value on Fisher's exact test.}
    \item{label}{A character string containing meta-data in sentence formatting,
    for convenient use in, e.g., plot captions.}
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
\code{\link{findAssociatedClones}}
\code{\link{buildAssociatedClusterNetwork}}
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

associated_seqs <- findAssociatedSeqs(
  file_list = file.path(data_dir, paste0(1:num_samples, ".rds")),
  input_type = "rds",
  group_ids = groups,
  min_seq_length = NULL,
  drop_matches = NULL,
  groups = c(0, 1),
  seq_col = "seq",
  outfile = NULL)
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
