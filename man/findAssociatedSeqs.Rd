\name{findAssociatedSeqs}
\alias{findAssociatedSeqs}
\alias{findAssociatedSeqs2}
\title{
    Identify TCR/BCR Sequences Associated With a Binary Variable
}
\description{
    Part of the workflow \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}.

    Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data and a binary variable of interest such as a disease condition, treatment or clinical outcome, identify receptor sequences that exhibit a statistically significant difference in frequency between the two levels of the binary variable.

    \code{findAssociatedSeqs} is designed for use when each sample is stored in a separate file. \code{findAssociatedSeqs2} is designed for use with a single data frame containing all samples.
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
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the data for a sample, with observations indexed by row. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{input_type}{
  A character string specifying the file format of the sample data files. Options are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"}, \code{"rds"} and \code{"rda"}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each sample's data frame within its respective Rdata file.  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sample_ids}{
  Deprecated. Does nothing.
}
  \item{subject_ids}{
  A character or numeric vector of subject IDs, whose length matches that of \code{file_list}. Only relevant when the binary variable of interest is subject-specific and multiple samples belong to the same subject.
}
  \item{group_ids}{
  A character or numeric vector of group IDs containing exactly two unique values and with length matching that of \code{file_list}. The two groups correspond to the two values of the binary variable of interest.
}
  \item{groups}{
  Deprecated. Does nothing.
}
  \item{seq_col}{
    Specifies the column of each sample's data frame containing the TCR/BCR sequences. Accepts either the column name or column index.
}
  \item{freq_col}{
    Optional. Specifies the column of each sample's data frame containing the clone frequency (i.e., clone count divided by the sum of the clone counts across all clones in the sample). Accepts either the column name or column index. If this argument is specified, the maximum clone frequency (across all samples) for each associated sequence will be included in the content of the \code{label} variable of the returned data frame.
}
  \item{min_seq_length}{
    Controls the minimum TCR/BCR sequence length considered when searching for associated sequences. Passed to \code{\link{filterInputData}}.
}
  \item{drop_matches}{
    Passed to \code{\link{filterInputData}}. Accepts a regular expression (see \code{\link[base]{regex}}). Checks TCR/BCR sequences for a pattern match using \code{\link{grep}}. Those returning a match are excluded from consideration as associated sequences. It is recommended to filter out sequences containing special characters that are invalid for use in file names. By default, sequences containing either of the characters \code{*} or \code{_} are dropped.
}
  \item{min_sample_membership}{
    Controls the minimum number of samples in which a TCR/BCR sequence must be present in order to be considered when searching for associated sequences. Setting this value to \code{NULL} bypasses the check.
}
  \item{pval_cutoff}{
    Controls the P-value cutoff below which an association is detected by Fisher's exact test (see details).
}
  \item{outfile}{
    A file path for saving the output (using \code{write.csv}).
}
\item{data}{
    A data frame containing the combined AIRR-seq data for all samples.
}
\item{sample_col}{
    The column of \code{data} containing the sample IDs. Accepts the column name or column index.
}
\item{subject_col}{
    Optional. The column of \code{data} containing the subject IDs. Accepts the column name or column index. Only relevant when the binary variable of interest is subject-specific and multiple samples belong to the same subject.
}
\item{group_col}{
    The column of \code{data} containing the group IDs. Accepts the column name or column index. The groups correspond to the two values of the binary variable of interest. Thus there should be exactly two unique values in this column.
}
}
\details{
    The TCR/BCR sequences from all samples are first filtered according to minimum sequence length and sequence content based on the specified values in \code{min_seq_length} and \code{drop_matches}, respectively. The sequences are further filtered based on sample membership, removing sequences appearing in fewer than \code{min_sample_membership} samples.

    For each remaining TCR/BCR sequence, a P-value is computed for Fisher's exact test of independence between the binary variable of interest and the presence of the sequence within a repertoire. The samples/subjects are divided into two groups based on the levels of the binary variable. If subject IDs are provided, then the test is based on the number of subjects in each group for whom the sequence appears in one of their samples. Without subject IDs, the test is based on the number of samples possessing the sequence in each group.

  Fisher's exact test is performed using \code{\link[stats]{fisher.test}}. TCR/BCR sequences with a \eqn{P}-value below \code{pval_cutoff} are sorted by \eqn{P}-value and returned along with some additional information.

    The returned ouput is intended for use with the \code{\link{findAssociatedClones}} function. See the \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters} package vignette for more details. The vignette can be accessed offline using \code{vignette("associated_clusters")}.
}
\value{
    A data frame containing the TCR/BCR sequences found to be associated with the binary variable using Fisher's exact test (see details). Each row corresponds to a unique TCR/BCR sequence and includes the following variables:
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
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters vignette}
}

\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{findAssociatedClones}}
\code{\link{buildAssociatedClusterNetwork}}
}
\examples{

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

# View first few rows of data for sample 1
head(readRDS(file.path(dir_input_samples, "Sample1.rds")))


## Step 1: Find Associated Sequences ##
# input files (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))

# group labels for the samples
group_labels <- c(rep("reference", samples_c), rep("comparison", samples_t))

# search for associated sequences using Fisher's exact test
associated_seqs <- findAssociatedSeqs(
  file_list = input_files, input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  min_seq_length = NULL, drop_matches = NULL,
  min_sample_membership = 0,
  pval_cutoff = 0.1, outfile = NULL)

# view first few rows of output
head(associated_seqs)


# clean up temp directory
unlink(list.dirs(data_dir, full.names = TRUE),
       recursive = TRUE)

}
