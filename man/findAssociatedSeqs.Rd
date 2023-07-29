\name{findAssociatedSeqs}
\alias{findAssociatedSeqs}
\alias{findAssociatedSeqs2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Identify Receptor Sequences Associated With a Binary Variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data and a sample-level or subject-level binary variable of interest such as a disease condition, treatment or clinical outcome, identify receptor sequences that exhibit a statistically significant difference in frequency between the two levels of the binary variable using Fisher's exact test.

    \code{findAssociatedSeqs} is designed for use when each sample is contained in a separate data frame, with each data frame stored in a separate file. \code{findAssociatedSeqs2} is designed for use with a single data frame containing all samples.
}
\usage{
findAssociatedSeqs(
  ## Input ##
  file_list,
  input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids = NULL,
  subject_ids = NULL,
  group_ids,
  groups = NULL,
  seq_col,
  freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7,
  drop_matches = "[*|_]",
  min_sample_membership = 5,
  pval_cutoff = 0.05,

  ## Output ##
  outfile = "associated_seqs.csv")


findAssociatedSeqs2(
  ## Input ##
  data,
  seq_col,
  sample_col,
  subject_col = sample_col,
  group_col,
  groups = NULL,
  freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7,
  drop_matches = "[*|_]",
  min_sample_membership = 5,
  pval_cutoff = 0.05,

  ## Ouptut ##
  outfile = "associated_seqs.csv")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{\link{loadDataFromFileList}} (one data frame per sample).
}
  \item{input_type}{
%%     ~~Describe \code{input_type} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
%%     ~~Describe \code{data_symbols} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sample_ids}{
%%     ~~Describe \code{sample_ids} here~~
  Deprecated. Does nothing.
}
  \item{subject_ids}{
%%     ~~Describe \code{subject_ids} here~~
  A character or numeric vector of subject IDs, whose length matches that of \code{file_list}.
  Only relevant when the binary variable of interest is subject-specific and multiple samples belong to the same subject.
}
  \item{group_ids}{
%%     ~~Describe \code{group_ids} here~~
  A character or numeric vector of group IDs, which contains exactly two unique
  values and whose length matches that of \code{file_list}. The two groups correspond to the two values of the binary variable of interest.
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
  Deprecated. Does nothing. This argument was previously used to specify the values of \code{group_ids} corresponding to the two values of the binary variable of interest. These values are now determined automatically from the unique values of \code{group_ids}.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The name or index of of the column containing the clone sequences in the input data for each sample.
}
  \item{freq_col}{
%%     ~~Describe \code{freq_col} here~~
    Optional. The name or index of the column containing the clone frequency/fraction (i.e., measure of clonal abundance normalized based on the total count in the sample) in the input data for each sample. If supplied, the maximum clone frequency (across all samples) for each associated sequence will be included in the content of the \code{label} variable of the returned data frame.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    Passed to \code{\link{filterInputData}}.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{\link{filterInputData}}.
}
  \item{min_sample_membership}{
%%     ~~Describe \code{min_sample_membership} here~~
    Only sequences that appear in at least this many samples will be considered. Setting this value to \code{NULL} bypasses the check.
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
    A data frame containing the combined AIRR-seq data for all samples,
    with variables indexed by column and observations indexed by row.
}
\item{sample_col}{
    The name or index of the column of \code{data} containing the sample IDs.
}
\item{subject_col}{
    Optional. The name or index of the column of \code{data} containing the subject IDs. Only relevant when the binary variable of interest is subject-specific and multiple samples belong to the same subject.
}
\item{group_col}{
    The name or index of the column of \code{data} containing the group IDs.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The list of unique receptor sequences is first filtered according to minimum sequence length and sequence content based on the specified values in \code{min_seq_length} and \code{drop_matches}, respectively. The list is then further filtered based on minimum sample membership, keeping
    only those receptor sequences that appear in at least \code{min_sample_membership}
    samples.

    For each remaining sequence, a P-value is computed for Fisher's exact test of independence between the binary variable of interest and the appearance of the sequence within a repertoire. The test is performed using
    \code{stats::fisher.test}, based on the number of samples possessing the
    sequence in each of the two groups. If subject IDs are provided, then the test is instead based on the number of subjects possessing the sequence in one or more samples. Only sequences
    with a \eqn{P}-value below \code{pval_cutoff} will be retained. The resulting
    sequences are sorted by \eqn{P}-value in increasing order and returned
    along with some basic meta-data.

    The returned ouput is intended for use with the \code{findAssociatedClones} function. See the Finding Associated Clusters package vignette for more details. The vignette can be accessed using \code{vignette("associated_clusters")}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the receptor sequences found to be associated with the binary variable based on Fisher's exact test using the specified P-value cutoff. Each row corresponds to a unique sequence and includes the following variables:
    \item{ReceptorSeq}{The unique receptor sequence}
    \item{fisher_pvalue}{The P-value on Fisher's exact test for independence between the receptor sequence and the binary variable of interest}
    \item{shared_by_n_samples}{The number of samples in which the sequence was observed}
    \item{samples_g0}{Of the samples in which the sequence was observed, the number of samples belonging to the first group}
    \item{samples_g1}{Of the samples in which the sequence was observed, the number of samples belonging to the second group}
    \item{shared_by_n_subjects}{The number of subjects in which the sequence was observed (only present if subject IDs are specified)}
    \item{subjects_g0}{Of the subjects in which the sequence was observed, the number of subjects belonging to the first group (only present if subject IDs are specified)}
    \item{subjects_g1}{Of the subjects in which the sequence was observed, the number of subjects belonging to the second group (only present if subject IDs are specified)}
    \item{label}{A character string summarizing the above information. Also includes the maximum in-sample clone frequency across all samples, if available}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. https://doi.org/10.3389/fimmu.2023.1181825

\url{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825}

\url{https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-}
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{findAssociatedClones}}
\code{\link{buildAssociatedClusterNetwork}}
\code{vignette("associated_clusters")}
}
\examples{

\dontrun{
## Generate some toy data for demonstration ##
# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_input_samples <- file.path(data_dir, "input_samples")
dir.create(dir_input_samples, showWarnings = FALSE)

# Number of samples by control/treatment group
samples_c <- samples_t <- 15
samples <- samples_c + samples_t
sample_size <- 30 # (seqs per sample)

# sequences (first five are chosen to be associated with treatment)
base_seqs <- c("CASSGAYEQYF", "CSVDLGKGNNEQFF",
               "CASSIEGQLSTDTQYF",
               "CASSEEGQLSTDTQYF",
               "CASSPEGQLSTDTQYF",
               "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF")

# relative generation probabilities by control/treatment group
pgen_c <- matrix(rep(c(rep(1, 5), rep(30, 3)), times = samples_c),
                 nrow = samples_c, byrow = TRUE)
pgen_t <- matrix(rep(c(1, 1, rep(1/3, 3), rep(2, 3)), times = samples_t),
                 nrow = samples_t, byrow = TRUE)
pgen <- rbind(pgen_c, pgen_t)

# Simulate the data
simulateToyData(
  samples = samples, sample_size = sample_size,
  prefix_length = 1, prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = base_seqs, affix_probs = pgen, num_edits = 0,
  output_dir = dir_input_samples, no_return = TRUE)


## Find Associated Sequences ##
# input files (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))

# group labels for the samples
group_labels <- c(rep("reference", samples_c), rep("comparison", samples_t))

associated_seqs <- findAssociatedSeqs(
  file_list = input_files, input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  min_seq_length = NULL, drop_matches = NULL,
  min_sample_membership = 0,
  pval_cutoff = 0.1,outfile = NULL)
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
