\name{simulateToyData}
\alias{simulateToyData}
\title{
Generate Toy AIRR-Seq Data
}
\description{
Generates toy data that can be used to test or demonstrate the behavior of
functions in the \code{NAIR} package.  Created as a lightweight tool for use in
tests, examples and vignettes. This function is not intended to simulate realistic
data.
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
  no_return = FALSE
)
}
\arguments{
  \item{samples}{
    The number of distinct samples to include in the data.
}
  \item{chains}{
    The number of chains (either 1 or 2) for which to generate receptor sequences.
}
  \item{sample_size}{
    The number of observations to generate per sample.
}
  \item{prefix_length}{
    The length of the random prefix generated for each observed sequence.
    Specifically, the number of elements of \code{prefix_chars} that are sampled
    with replacement and concatenated to form each prefix.
}
  \item{prefix_chars}{
    A character vector containing characters or strings from which to sample when
    generating the prefix for each observed sequence.
}
  \item{prefix_probs}{
    A numeric matrix whose column dimension matches the length of \code{prefix_chars}
    and with row dimension matching the value of \code{samples}. The \eqn{i}th row
    specifies the relative probability weights assigned to each element of
    \code{prefix_chars} when sampling to form the prefix for each sequence in the
    \eqn{i}th sample.
}
  \item{affixes}{
    A character vector containing characters or strings from which to sample when
    generating the suffix for each observed sequence.
}
  \item{affix_probs}{
    A numeric matrix whose column dimension matches the length of \code{affixes}
    and with row dimension matching the value of \code{samples}. The \eqn{i}th
    row specifies the relative probability weights assigned to each element of
    \code{affixes} when sampling to form the suffix for each sequence in the
    \eqn{i}th sample.
}
  \item{num_edits}{
    A nonnegative integer specifying the number of random edit operations to
    perform on each observed sequence after its initial generation.
}
  \item{edit_pos_probs}{
    A function that accepts a nonnegative integer (the character length of a
    sequence) as its argument and returns a vector of this length containing
    probability weights. Each time an edit operation is performed on a sequence,
    the character position at which to perform the operation is randomly
    determined according to the probabilities given by this function.
}
  \item{edit_ops}{
    A character vector specifying the possible operations that can be performed
    for each edit. The default value includes all valid operations (insertion,
    deletion, transmutation).
}
  \item{edit_probs}{
    A numeric vector of the same length as \code{edit_ops}, specifying the relative
    probability weights assigned to each edit operation.
}
  \item{new_chars}{
    A character vector containing characters or strings from which to sample when
    performing an insertion edit operation.
}
  \item{new_probs}{
    A numeric matrix whose column dimension matches the length of \code{new_chars}
    and with row dimension matching the value of \code{samples}. The \eqn{i}th row
    specifies, for the \eqn{i}th sample, the relative probability weights assigned
    to each element of \code{new_chars} when performing a transmutation or insertion
    as a random edit operation.
}
  \item{output_dir}{
    An optional character string specifying a file directory to save the generated
    data. One file will be generated per sample.
}
  \item{no_return}{
    A logical flag that can be used to prevent the function from returning the
    generated data. If \code{TRUE}, the function will instead return \code{TRUE}
    once all processes are complete.
}
}
\details{
    Each observed sequence is obtained by separately generating a prefix and
    suffix according to the specified settings, then joining the two and performing
    sequential rounds of edit operations randomized according to the user's
    specifications.

    Count data is generated for each observation; note that this count data is
    generated independently from the observed sequences and has no relationship
    to them.
}
\value{
    If \code{no_return = FALSE} (the default), a \code{data.frame} whose contents depend
    on the value of the \code{chains} argument.

    For \code{chains = 1}, the data frame contains the following variables:
    \item{CloneSeq}{The "receptor sequence" for each observation.}
    \item{CloneFrequency}{The "clone frequency" for each observation (clone count
    as a proportion of the aggregate clone count within each sample).}
    \item{CloneCount}{The "clone count" for each observation.}
    \item{SampleID}{The sample ID for each observation.}

    For \code{chains = 2}, the data frame contains the following variables:
    \item{AlphaSeq}{The "alpha chain" receptor sequence for each observation.}
    \item{AlphaSeq}{The "beta chain" receptor sequence for each observation.}
    \item{UMIs}{The "unique molecular identifier count" for each observation.}
    \item{Count}{The "count" for each observation.}
    \item{SampleID}{The sample ID for each observation.}

    If \code{no_return = TRUE}, the function returns \code{TRUE} upon completion.
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
set.seed(42)

# Bulk data from two samples
dat1 <- simulateToyData()

# Single-cell data with alpha and beta chain sequences
dat2 <- simulateToyData(chains = 2)

# Write data to file, return nothing
simulateToyData(sample_size = 500,
                num_edits = 10,
                no_return = TRUE,
                output_dir = tempdir())

# Example customization
dat4 <-
  simulateToyData(
    samples = 5,
    sample_size = 50,
    prefix_length = 0,
    prefix_chars = "",
    prefix_probs = matrix(1, nrow = 5),
    affixes = c("CASSLGYEQYF", "CASSLGETQYF",
                "CASSLGTDTQYF", "CASSLGTEAFF",
                "CASSLGGTEAFF", "CAGLGGRDQETQYF",
                "CASSQETQYF", "CASSLTDTQYF",
                "CANYGYTF", "CANTGELFF",
                "CSANYGYTF"),
    affix_probs = matrix(1, ncol = 11, nrow = 5),
  )

## Simulate 30 samples with a mix of public/private sequences ##
samples <- 30
sample_size <- 30 # (seqs per sample)
base_seqs <- c(
  "CASSIEGQLSTDTQYF", "CASSEEGQLSTDTQYF", "CASSSVETQYF",
  "CASSPEGQLSTDTQYF", "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF",
  "CASSLTSGYNEQFF", "CASSETGYNEQFF", "CASSLTGGNEQFF", "CASSYLTGYNEQFF",
  "CASSLTGNEQFF", "CASSLNGYNEQFF", "CASSFPWDGYGYTF", "CASTLARQGGELFF",
  "CASTLSRQGGELFF", "CSVELLPTGPLETSYNEQFF", "CSVELLPTGPSETSYNEQFF",
  "CVELLPTGPSETSYNEQFF", "CASLAGGRTQETQYF", "CASRLAGGRTQETQYF",
  "CASSLAGGRTETQYF", "CASSLAGGRTQETQYF", "CASSRLAGGRTQETQYF",
  "CASQYGGGNQPQHF", "CASSLGGGNQPQHF", "CASSNGGGNQPQHF", "CASSYGGGGNQPQHF",
  "CASSYGGGQPQHF", "CASSYKGGNQPQHF", "CASSYTGGGNQPQHF",
  "CAWSSQETQYF", "CASSSPETQYF", "CASSGAYEQYF", "CSVDLGKGNNEQFF")
# Relative generation probabilities
pgen <- cbind(
  stats::toeplitz(0.6^(0:(sample_size - 1))),
  matrix(1, nrow = samples, ncol = length(base_seqs) - samples))
dat5 <-
  simulateToyData(
    samples = samples,
    sample_size = sample_size,
    prefix_length = 1,
    prefix_chars = c("", ""),
    prefix_probs = cbind(rep(1, samples), rep(0, samples)),
    affixes = base_seqs,
    affix_probs = pgen,
    num_edits = 0
  )

## Simulate 30 samples from two groups (treatment/control) ##
samples_c <- samples_t <- 15 # Number of samples by control/treatment group
samples <- samples_c + samples_t
sample_size <- 30 # (seqs per sample)
base_seqs <- # first five are associated with treatment
  c("CASSGAYEQYF", "CSVDLGKGNNEQFF", "CASSIEGQLSTDTQYF",
    "CASSEEGQLSTDTQYF", "CASSPEGQLSTDTQYF",
    "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF")
# Relative generation probabilities by control/treatment group
pgen_c <- matrix(rep(c(rep(1, 5), rep(30, 3)), times = samples_c),
                 nrow = samples_c, byrow = TRUE)
pgen_t <- matrix(rep(c(1, 1, rep(1/3, 3), rep(2, 3)), times = samples_t),
                 nrow = samples_t, byrow = TRUE)
pgen <- rbind(pgen_c, pgen_t)
dat6 <-
  simulateToyData(
    samples = samples,
    sample_size = sample_size,
    prefix_length = 1,
    prefix_chars = c("", ""),
    prefix_probs =
      cbind(rep(1, samples), rep(0, samples)),
    affixes = base_seqs,
    affix_probs = pgen,
    num_edits = 0
  )


\dontshow{
# clean up temp directory
file.remove(
  file.path(tempdir(), paste0("Sample", 1:2, ".rds"))
)
}
}
