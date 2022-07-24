\name{getSimilarClones}
\alias{getSimilarClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Find Clones Similar to a Target by Receptor Sequence
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given immune repertoire sequence data and a target receptor sequence, return
    the subset of the data consiting of those clones whose receptor sequences are
    within a specified distance of the target.
}
\usage{
getSimilarClones(target_seq, data, clone_col, sample_col = NULL,
                 dist_type = "hamming", max_dist = 2, drop_chars = "[*|_]")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target_seq}{
%%     ~~Describe \code{target_seq} here~~
    A character string containing the target receptor sequence.
}
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones
      indexed by row.
}
  \item{clone_col}{
%%     ~~Describe \code{clone_col} here~~
      The column name or number of `data` containing the receptor sequences (e.g.,
      TCR CDR3 nucleotide sequence or amino acid sequence).
}
  \item{sample_col}{
%%     ~~Describe \code{sample_col} here~~
      (Optional) The column name or number of `data` containing the sample ID.
      If supplied, only clones from samples that contain the target sequence
      will be included. The sample IDs can be either character or numeric.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
      The distance metric to use as a measure of similarity between two receptor
      sequences. Valid options are `"hamming"` (the default) and `"levenshtein"`
      (i.e., edit distance).
}
  \item{max_dist}{
%%     ~~Describe \code{max_dist} here~~
      A numeric value used as a cutoff; clones whose receptor sequences have
      distance from the target sequence exceeding `max_dist` (where the distance
      type is specified by `dist_type`) will be removed.
}
  \item{drop_chars}{
      A regular expression; clones whose receptor sequences match the specified
      expression will be removed. By default, this checks for sequences
      containing any of the special characters `*`, `_`, or `|`.
%%     ~~Describe \code{drop_chars} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the same columns as `data` and those rows of `data`
    for which the sequence in `clone_col` has `dist_type` distance from
    `target_seq` of no greater than `max_dist` (and does not match the expression
    specified by `drop_chars`). If `sample_col` is non-null, only clones belonging
    to samples that possess the target sequence will be included.
}
\references{
%% ~put references to the literature/web site here ~
    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Create some data
data <- data.frame(
  clone_seq = c("ATCG", "TTCG", "TCG", rep("ATCC", 2), "ATGG", "TGG"),
  clone_count = rep(1, 7),
  clone_freq = rep(1 / 7, 7),
  # group/label variable 1
  time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
  # group/label variable 2
  subject_id = c(rep(1, 5), rep(2, 2))
)
target <- "ATCG"

# Default usage (Hamming distance, max distance 1)
data_target <-
  getSimilarClones(target, data, "clone_seq")

# Max Hamming distance of 2
# No sequences have a Hamming distance of 2 from the target,
# so the output is the same as in the previous example
data_target2 <-
  getSimilarClones(target, data, "clone_seq",
                   max_dist = 2)

# Levenshtein (edit) distance
data_target_lev <-
  getSimilarClones(target, data, "clone_seq",
                   dist_type = "levenshtein")

# Levenshtein (edit) distance, max distance 2
data_target_lev2 <-
  getSimilarClones(target, data, "clone_seq",
                   dist_type = "levenshtein",
                   max_dist = 2)

# Keeping only the subjects with the target sequence
data_subjects_lev2 <-
  getSimilarClones(target, data, "clone_seq",
                   sample_col = "subject_id",
                   dist_type = "levenshtein",
                   max_dist = 2)

# Keeping only the time points with the target sequence
data_timepoints_lev2 <-
  getSimilarClones(target, data, "clone_seq",
                   sample_col = "time_point",
                   dist_type = "levenshtein",
                   max_dist = 2)
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
