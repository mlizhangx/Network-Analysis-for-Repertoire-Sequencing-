\name{findAssociatedClones}
\alias{findAssociatedClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Gather Clones in Public Neighborhoods Around Associated Sequences
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple samples of RepSeq data from two groups and a list of receptor
    sequences associated to the comparison group (e.g., obtained using
    \code{findAssociatedSeqs}), this function does the following:
    for each associated sequence,
    from all samples that possess the associated sequence, gather all clones
    whose receptor sequences lie in a neighborhood of the associated sequence.

    The public neighborhood for each associated sequence is stored in its own data frame
    and saved to its own file, in preparation for use with
    \code{buildAssociatedClusterNetwork}.
}
\usage{
findAssociatedClones(

  ## Input ##
  file_list, input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids = 1:length(file_list),
  subject_ids = sample_ids,
  group_ids,
  seq_col,

  ## Search Criteria ##
  assoc_seqs,
  nbd_radius = 1,
  dist_type = "hamming",
  min_seq_length = 6,
  drop_matches = "[*|_]",

  ## Output ##
  subset_cols = NULL,
  output_dir =
    file.path(getwd(), "associated_neighborhoods"),
  output_type = "csv"
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} (one file per sample).
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
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of each sample's data frame that contains the receptor sequences to associate with the comparison group.
}
  \item{assoc_seqs}{
%%     ~~Describe \code{assoc_seqs} here~~
    A character vector containing the receptor sequences associated to the comparison group.
}
  \item{nbd_radius}{
%%     ~~Describe \code{nbd_radius} here~~
    The maximum distance between an element of \code{assoc_seqs} and other clones within its neighborhood.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    The distance metric used to measure the distance between sequences in \code{nbd_radius}. Valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    Passed to \code{filterInputData()} when loading each sample.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{filterInputData()}.
}
  \item{subset_cols}{
%%     ~~Describe \code{subset_cols} here~~
    Passed to \code{filterInputData()}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
  The output directory; a valid output directory is required, since results are not returned in R.
}
  \item{output_type}{
%%     ~~Describe \code{output_dir} here~~
  A character string specifying the file format to use when writing each file. Valid options include \code{"csv"}, \code{"tsv"}, \code{"rds"} and \code{"rda"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
   Each neighborhood's RepSeq data file is saved using the corresponding associated sequence
   for the filename, with the appropriate file extension appended based on the
   value of \code{output_type} (e.g., the neighborhood data for the sequence
   \code{"CASSGAYEQYF"} would be saved as \code{"CASSGAYEQYF.csv"} if
   \code{output_type = "csv"}).
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
