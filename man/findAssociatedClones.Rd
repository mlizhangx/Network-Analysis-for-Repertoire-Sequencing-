\name{findAssociatedClones}
\alias{findAssociatedClones}
\title{
    Identify TCR/BCR Clones in a Neighborhood Around Each Associated Sequence
}
\description{
Part of the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}.
Intended for use following \code{\link[=findAssociatedSeqs]{findAssociatedSeqs()}}
and prior to
\code{\link[=buildAssociatedClusterNetwork]{buildAssociatedClusterNetwork()}}.

Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing
(AIRR-Seq) data and a vector of associated sequences, identifies for each
associated sequence a global "neighborhood" comprised of clones with TCR/BCR
sequences similar to the associated sequence.
}
\usage{
findAssociatedClones(

  ## Input ##
  file_list, input_type,
  data_symbols = NULL,
  header, sep, read.args,
  sample_ids =
    paste0("Sample", 1:length(file_list)),
  subject_ids = NULL,
  group_ids,
  seq_col,
  assoc_seqs,

  ## Neighborhood Criteria ##
  nbd_radius = 1,
  dist_type = "hamming",
  min_seq_length = 6,
  drop_matches = NULL,

  ## Output ##
  subset_cols = NULL,
  output_dir,
  output_type = "rds",
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
  A character string specifying the file format of the sample data files. Options
  are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"}, \code{"rds"} and
  \code{"rda"}.
  Passed to \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each sample's data
  frame within its respective Rdata file.
  Passed to \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
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
  A character or numeric vector of sample IDs, whose length matches that of
  \code{file_list}. Each entry is assigned as the sample ID to the corresponding
  entry of \code{file_list}.
}
  \item{subject_ids}{
  An optional character or numeric vector of subject IDs, whose length matches
  that of \code{file_list}. Used to assign a subject ID to each sample.
}
  \item{group_ids}{
  A character or numeric vector of group IDs whose length matches that of
  \code{file_list}. Used to assign each sample to a group. The two groups
  represent the levels of the binary variable of interest.
}
  \item{seq_col}{
    Specifies the column of each sample's data frame containing the TCR/BCR
    sequences.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
  \item{assoc_seqs}{
    A character vector containing the TCR/BCR sequences associated with the
    binary variable of interest.
}
  \item{nbd_radius}{
    The maximum distance (based on \code{dist_type}) between an associated sequence
    and other TCR/BCR sequences belonging to its neighborhood. Lower values require
    sequences to be more similar to an associated sequence in order to belong to
    its neighborhood.
}
  \item{dist_type}{
    Specifies the function used to quantify the similarity between sequences. The
    similarity between two sequences determines their pairwise distance, with
    greater similarity corresponding to shorter distance. Valid options are
    \code{"hamming"} (the default), which
    uses \code{\link[=hamDistBounded]{hamDistBounded()}},
    and \code{"levenshtein"}, which uses
    \code{\link[=levDistBounded]{levDistBounded()}}.
}
  \item{min_seq_length}{
    Clones with TCR/BCR sequences below this length will be removed. Passed to
    \code{\link[=filterInputData]{filterInputData()}} when loading each sample.
}
  \item{drop_matches}{
    Passed to \code{\link[=filterInputData]{filterInputData()}}.
    Accepts a character string containing a regular expression
    (see \code{\link[base]{regex}}).
    Checks TCR/BCR sequences for a pattern match using \code{\link[=grep]{grep()}}.
    Those returning a match are dropped.
    By default, sequences containing any of the characters \code{*}, \code{|}
    or \code{_} are dropped.
}
  \item{subset_cols}{
    Controls which columns of the AIRR-Seq data from each sample are included in
    the output.
    Accepts a character vector of column names
    or a numeric vector of column indices.
    The default \code{NULL} includes all columns. Passed to
    \code{\link[=filterInputData]{filterInputData()}}.
}
  \item{output_dir}{
  A file path to a directory for saving the output. A valid output directory is
  required, since no output is returned in R. The specified directory will be
  created if it does not already exist.
}
  \item{output_type}{
  A character string specifying the file format to use for saving the output.
  Valid options are
  \code{"rda"}, \code{"csv"}, \code{"csv2"}, \code{"tsv"} and\code{"table"}.
  For \code{"rda"}, data frames are named \code{data} in the R environment.
  For the remaining options, \code{\link[utils:write.table]{write.table()}} is
  called with \code{row.names = TRUE}.
}
  \item{verbose}{
  Logical. If \code{TRUE}, generates messages about the tasks
  performed and their progress, as well as relevant properties of intermediate
  outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
}
\details{
For each associated sequence, its neighborhood is defined to include all
clones with TCR/BCR sequences that are sufficiently similar to the associated
sequence. The arguments \code{dist_type} and \code{nbd_radius} control how the
similarity is measured and the degree of similarity required for neighborhood
membership.

For each associated sequence, a data frame is saved to an individual file.
The data frame contains one row for each clone in the associated sequence's
neighborhood (from all samples). It includes variables for sample ID, group ID
and (if provided) subject ID, as well as variables from the AIRR-Seq data.

The files saved by this function are intended for use with
\code{\link[=buildAssociatedClusterNetwork]{buildAssociatedClusterNetwork()}}. See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}
article on the package website for more details.
}
\value{
Returns \code{TRUE}, invisibly. The function is called for its side effects.
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
\code{\link[=findAssociatedSeqs]{findAssociatedSeqs()}}
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

