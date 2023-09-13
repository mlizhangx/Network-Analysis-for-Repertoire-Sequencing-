\name{findAssociatedSeqs}
\alias{findAssociatedSeqs}
\alias{findAssociatedSeqs2}
\title{
    Identify TCR/BCR Sequences Associated With a Binary Variable
}
\description{
Part of the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}.

Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing
(AIRR-Seq) data and a binary variable of interest such as a disease condition,
treatment or clinical outcome, identify receptor sequences that exhibit a
statistically significant difference in frequency between the two levels of the
binary variable.

\code{findAssociatedSeqs()} is designed for use when each sample is stored in a
separate file. \code{findAssociatedSeqs2()} is designed for use with a single data
frame containing all samples.
}
\usage{
findAssociatedSeqs(
  ## Input ##
  file_list,
  input_type,
  data_symbols = NULL,
  header, sep, read.args,
  sample_ids = deprecated(),
  subject_ids = NULL,
  group_ids,
  groups = deprecated(),
  seq_col,
  freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7,
  drop_matches = "[*|_]",
  min_sample_membership = 5,
  pval_cutoff = 0.05,

  ## Output ##
  outfile = NULL,
  verbose = FALSE
)


findAssociatedSeqs2(
  ## Input ##
  data,
  seq_col,
  sample_col,
  subject_col = sample_col,
  group_col,
  groups = deprecated(),
  freq_col = NULL,

  ## Search Criteria ##
  min_seq_length = 7,
  drop_matches = "[*|_]",
  min_sample_membership = 5,
  pval_cutoff = 0.05,

  ## Ouptut ##
  outfile = NULL,
  verbose = FALSE
)

}
\arguments{
  \item{file_list}{
  A character vector of file paths, or a list containing
  \code{\link[=connection]{connections}} and file paths.
  Each element corresponds to a single file containing the data
  for a single sample.
  Passed to \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{input_type}{
  A character string specifying the file format of the sample data files.
  Options are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"},
  \code{"rds"} and \code{"rda"}.
  Passed to \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each sample's data
  frame within its respective Rdata file.  Passed to
  \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{header}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify a non-default value of the \code{header}
    argument to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
}
  \item{sep}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify a non-default value of the \code{sep}
    argument to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
}
  \item{read.args}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify non-default values of optional
    arguments to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
    Accepts a named list of argument values.
    Values of \code{header} and \code{sep}
    in this list take precedence over values specified via the \code{header}
    and \code{sep} arguments.
}
  \item{sample_ids}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
  Does nothing.
}
  \item{subject_ids}{
  A character or numeric vector of subject IDs, whose length matches that of
  \code{file_list}. Only relevant when the binary variable of interest is
  subject-specific and multiple samples belong to the same subject.
}
  \item{group_ids}{
  A character or numeric vector of group IDs containing exactly two unique values
  and with length matching that of \code{file_list}. The two groups correspond to
  the two values of the binary variable of interest.
}
  \item{groups}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
  Does nothing.
}
  \item{seq_col}{
    Specifies the column of each sample's data frame containing the TCR/BCR sequences.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
  \item{freq_col}{
    Optional. Specifies the column of each sample's data frame containing the clone
    frequency (i.e., clone count divided by the sum of the clone counts across all
    clones in the sample).
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
    If this
    argument is specified, the maximum clone frequency (across all samples) for
    each associated sequence will be included in the content of the \code{label}
    variable of the returned data frame.
}
  \item{min_seq_length}{
    Controls the minimum TCR/BCR sequence length considered when searching for
    associated sequences. Passed to \code{\link[=filterInputData]{filterInputData()}}.
}
  \item{drop_matches}{
    Passed to \code{\link[=filterInputData]{filterInputData()}}.
    Accepts a character string containing
    a regular expression (see \code{\link[base]{regex}}). Checks TCR/BCR sequences
    for a pattern match using \code{\link[=grep]{grep()}}. Those returning a match are
    excluded from consideration as associated sequences. It is recommended to
    filter out sequences containing special characters that are invalid for use
    in file names. By default, sequences containing any of the characters
    \code{*}, \code{|} or \code{_} are dropped.
}
  \item{min_sample_membership}{
    Controls the minimum number of samples in which a TCR/BCR sequence must be
    present in order to be considered when searching for associated sequences.
    Setting this value to \code{NULL} bypasses the check.
}
  \item{pval_cutoff}{
    Controls the P-value cutoff below which an association is detected by Fisher's
    exact test (see details).
}
  \item{outfile}{
    A file path for saving the output (using \code{\link[=write.csv]{write.csv()}}).
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
\item{data}{
    A data frame containing the combined AIRR-seq data for all samples.
}
\item{sample_col}{
    The column of \code{data} containing the sample IDs.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
\item{subject_col}{
    Optional. The column of \code{data} containing the subject IDs.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
    Only relevant when the binary variable of interest
    is subject-specific and multiple samples belong to the same subject.
}
\item{group_col}{
    The column of \code{data} containing the group IDs.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
      The groups correspond to the two values of the binary variable
    of interest. Thus there should be exactly two unique values in this column.
}
}
\details{
The TCR/BCR sequences from all samples are first filtered according to minimum
sequence length and sequence content based on the specified values in
\code{min_seq_length} and \code{drop_matches}, respectively. The sequences
are further filtered based on sample membership, removing sequences appearing
in fewer than \code{min_sample_membership} samples.

For each remaining TCR/BCR sequence, a P-value is computed for Fisher's exact
test of independence between the binary variable of interest and the presence
of the sequence within a repertoire. The samples/subjects are divided into two
groups based on the levels of the binary variable. If subject IDs are provided,
then the test is based on the number of subjects in each group for whom the
sequence appears in one of their samples. Without subject IDs, the test is
based on the number of samples possessing the sequence in each group.

Fisher's exact test is performed using
\code{\link[stats:fisher.test]{fisher.test()}}. TCR/BCR
sequences with a \eqn{P}-value below \code{pval_cutoff} are sorted by \eqn{P}-value
and returned along with some additional information.

The returned ouput is intended for use with the
\code{\link[=findAssociatedClones]{findAssociatedClones()}}
function. See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}
article on the package website.
}
\value{
    A data frame containing the TCR/BCR sequences found to be associated with
    the binary variable using Fisher's exact test (see details). Each row corresponds
    to a unique TCR/BCR sequence and includes the following variables:
    \item{ReceptorSeq}{
    The unique receptor sequence.
    }
    \item{fisher_pvalue}{
    The P-value on Fisher's exact test for independence between the receptor
    sequence and the binary variable of interest.
    }
    \item{shared_by_n_samples}{
    The number of samples in which the sequence was observed.
    }
    \item{samples_g0}{
    Of the samples in which the sequence was observed, the number of samples
    belonging to the first group.
    }
    \item{samples_g1}{
    Of the samples in which the sequence was observed, the number of samples
    belonging to the second group.
    }
    \item{shared_by_n_subjects}{
    The number of subjects in which the sequence was observed (only present if
    subject IDs are specified).
    }
    \item{subjects_g0}{
    Of the subjects in which the sequence was observed, the number of subjects
    belonging to the first group (only present if subject IDs are specified).
    }
    \item{subjects_g1}{
    Of the subjects in which the sequence was observed, the number of subjects
    belonging to the second group (only present if subject IDs are specified).
    }
    \item{max_freq}{
    The maximum clone frequency across all samples.
    Only present if \code{freq_col} is non-null.
    }
    \item{label}{
    A character string summarizing the above information. Also includes the
    maximum in-sample clone frequency across all samples, if available.
    }
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters article on package website}
}

\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link[=findAssociatedClones]{findAssociatedClones()}}
\code{\link[=buildAssociatedClusterNetwork]{buildAssociatedClusterNetwork()}}
}
\examples{
set.seed(42)

## Simulate 30 samples from two groups (treatment/control) ##
n_control <- n_treatment <- 15
n_samples <- n_control + n_treatment
sample_size <- 30 # (seqs per sample)
base_seqs <- # first five are associated with treatment
  c("CASSGAYEQYF", "CSVDLGKGNNEQFF", "CASSIEGQLSTDTQYF",
    "CASSEEGQLSTDTQYF", "CASSPEGQLSTDTQYF",
    "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF")
# Relative generation probabilities by control/treatment group
pgen_c <- matrix(rep(c(rep(1, 5), rep(30, 3)), times = n_control),
                 nrow = n_control, byrow = TRUE)
pgen_t <- matrix(rep(c(1, 1, rep(1/3, 3), rep(2, 3)), times = n_treatment),
                 nrow = n_treatment, byrow = TRUE)
pgen <- rbind(pgen_c, pgen_t)
simulateToyData(
  samples = n_samples,
  sample_size = sample_size,
  prefix_length = 1,
  prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, n_samples), rep(0, n_samples)),
  affixes = base_seqs,
  affix_probs = pgen,
  num_edits = 0,
  output_dir = tempdir(),
  no_return = TRUE
)

## Step 1: Find Associated Sequences ##
sample_files <-
  file.path(tempdir(),
            paste0("Sample", 1:n_samples, ".rds")
  )
group_labels <- c(rep("reference", n_control),
                  rep("comparison", n_treatment))
associated_seqs <-
  findAssociatedSeqs(
    file_list = sample_files,
    input_type = "rds",
    group_ids = group_labels,
    seq_col = "CloneSeq",
    min_seq_length = NULL,
    drop_matches = NULL,
    min_sample_membership = 0,
    pval_cutoff = 0.1
  )
head(associated_seqs[, 1:5])

## Step 2: Find Associated Clones ##
dir_step2 <- tempfile()
findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir_step2
)

## Step 3: Global Network of Associated Clusters ##
associated_clusters <-
  buildAssociatedClusterNetwork(
    file_list = list.files(dir_step2,
                           full.names = TRUE
    ),
    seq_col = "CloneSeq",
    size_nodes_by = 1.5,
    print_plots = TRUE
  )

\dontshow{
# clean up temp directory
file.remove(
  file.path(tempdir(),
            paste0("Sample", 1:n_samples, ".rds")
  )
)
unlink(dir_step2, recursive = TRUE)
}
}