\name{simulateToyData}
\alias{simulateToyData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate Toy RepSeq Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates toy data that can be used to test or demonstrate the behavior of functions in the \code{NAIR} package.
Created as a lightweight tool for use in tests, examples and vignettes.
This function is not intended to simulate realistic data.
}
\usage{
simulateToyData(
    samples = 2,
    chains = 1,
    sample_size = 100,
    prefix_length = 7,
    prefix_chars = c("G", "A", "T", "C"),
    prefix_probs = rbind(
      "sample1" = c(12, 4, 1, 1),
      "sample2" = c(4, 12, 1, 1)),
    affixes = c("AATTGG", "AATCGG", "AATTCG",
                "AATTGC", "AATTG", "AATTC"),
    affix_probs = rbind(
      "sample1" = c(10, 4, 2, 2, 1, 1),
      "sample2" = c(1, 1, 1, 2, 2.5, 2.5)),
    num_edits = 0,
    edit_pos_probs = function(seq_length) {
      stats::dnorm(seq(-4, 4, length.out = seq_length))
    },
    edit_ops = c("insertion", "deletion", "transmutation"),
    edit_probs = c(5, 1, 4),
    new_chars = prefix_chars,
    new_probs = prefix_probs,
    output_dir = NULL,
    no_return = FALSE,
    seed_value = 42
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{
%%     ~~Describe \code{samples} here~~
    The number of distinct samples to include in the data.
}
  \item{chains}{
%%     ~~Describe \code{chains} here~~
    The number of chains (either 1 or 2) for which to generate receptor sequences.
}
  \item{sample_size}{
%%     ~~Describe \code{sample_size} here~~
    The number of observations to generate per sample.
}
  \item{prefix_length}{
%%     ~~Describe \code{prefix_length} here~~
    The length of the random prefix generated for each observed sequence. Specifically, the number of elements of \code{prefix_chars} that are sampled with replacement and concatenated to form each prefix.
}
  \item{prefix_chars}{
%%     ~~Describe \code{prefix_chars} here~~
    A character vector containing characters or strings from which to sample when generating the prefix for each observed sequence.
}
  \item{prefix_probs}{
%%     ~~Describe \code{prefix_probs} here~~
    A numeric matrix whose column dimension matches the length of \code{prefix_chars} and with row dimension matching the value of \code{samples}. The \eqn{i}th row specifies the relative probability weights assigned to each element of \code{prefix_chars} when sampling to form the prefix for each sequence in the \eqn{i}th sample.
}
  \item{affixes}{
%%     ~~Describe \code{affixes} here~~
    A character vector containing characters or strings from which to sample when generating the suffix for each observed sequence.
}
  \item{affix_probs}{
%%     ~~Describe \code{affix_probs} here~~
    A numeric matrix whose column dimension matches the length of \code{affixes} and with row dimension matching the value of \code{samples}. The \eqn{i}th row specifies the relative probability weights assigned to each element of \code{affixes} when sampling to form the suffix for each sequence in the \eqn{i}th sample.
}
  \item{num_edits}{
%%     ~~Describe \code{num_edits} here~~
    A nonnegative integer specifying the number of random edit operations to perform on each observed sequence after its initial generation.
}
  \item{edit_pos_probs}{
%%     ~~Describe \code{edit_pos_probs} here~~
    A function that accepts a nonnegative integer (the character length of a sequence) as its argument and returns a vector of this length containing probability weights. Each time an edit operation is performed on a sequence, the character position at which to perform the operation is randomly determined according to the probabilities given by this function.
}
  \item{edit_ops}{
%%     ~~Describe \code{edit_ops} here~~
    A character vector specifying the possible operations that can be performed for each edit. The default value includes all valid operations (insertion, deletion, transmutation).
}
  \item{edit_probs}{
%%     ~~Describe \code{edit_probs} here~~
    A numeric vector of the same length as \code{edit_ops}, specifying the relative probability weights assigned to each edit operation.
}
  \item{new_chars}{
%%     ~~Describe \code{new_chars} here~~
    A character vector containing characters or strings from which to sample when performing an insertion edit operation.
}
  \item{new_probs}{
%%     ~~Describe \code{new_probs} here~~
    A numeric vector of the same length as \code{new_chars}, specifying the relative probability weights assigned to each possible character or string to be inserted.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    An optional character string specifying a file directory to save the generated data. One file will be generated per sample.
}
  \item{no_return}{
%%     ~~Describe \code{no_return} here~~
    A logical flag that can be used to prevent the function from returning the generated data. If \code{TRUE}, the function will instead return \code{TRUE} once all processes are complete.
}
  \item{seed_value}{
%%     ~~Describe \code{seed_value} here~~
    Passed to \code{set_seed} at the start of function execution.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    Each observed sequence is obtained by separately generating a prefix and suffix
    according to the specified settings, then joining the two and performing
    sequential rounds of edit operations randomized according to the user's
    specifications.

    Count data is generated for each observation; note that this count data is
    generated independently from the observed sequences and has no relationship
    to them.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If \code{no_return = FALSE} (the default), A data frame whose contents depend on the value of the \code{chains} argument.

    For \code{chains = 1}, the data frame contains the following variables:
    \item{CloneSeq}{The "receptor sequence" for each observation.}
    \item{CloneFrequency}{The "clone frequency" for each observation (clone count as a proportion of the aggregate clone count within each sample).}
    \item{CloneCount}{The "clone count" for each observation.}
    \item{SampleID}{The sample ID for each observation.}

    For \code{chains = 2}, the data frame contains the following variables:
    \item{AlphaSeq}{The "alpha chain" receptor sequence for each observation.}
    \item{AlphaSeq}{The "beta chain" receptor sequence for each observation.}
    \item{UMIs}{A "unique count" for each observation.}
    \item{Count}{The "count" for each observation.}
    \item{SampleID}{The sample ID for each observation.}
}
    If \code{no_return = FALSE}, the function returns \code{TRUE} upon completion.
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
    simulateToyData()
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
