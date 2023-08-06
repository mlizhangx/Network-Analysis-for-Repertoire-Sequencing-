\name{findAssociatedClones}
\alias{findAssociatedClones}
\title{
    Identify TCR/BCR Clones in a Neighborhood Around Each Associated Sequence
}
\description{
    Part of the workflow \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}. Intended for use following \code{\link{findAssociatedSeqs}} and prior to \code{\link{buildAssociatedClusterNetwork}}.

    Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data and a vector of associated sequences, identifies for each associated sequence a global "neighborhood" comprised of clones with TCR/BCR sequences similar to the associated sequence.
}
\usage{
findAssociatedClones(

  ## Input ##
  file_list, input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
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
  drop_matches = "[*|_]",

  ## Output ##
  subset_cols = NULL,
  output_dir =
    file.path(getwd(), "associated_neighborhoods"),
  output_type = "csv",
  verbose = FALSE
)
}
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the AIRR-Seq data for a sample, with observations indexed by row. Passed to \code{\link{loadDataFromFileList}}.
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
  A character or numeric vector of sample IDs, whose length matches that of \code{file_list}.
}
  \item{subject_ids}{
  An optional character or numeric vector of subject IDs, whose length matches that of \code{file_list}. Used to assign a subject ID to each sample.
}
  \item{group_ids}{
  A character or numeric vector of group IDs whose length matches that of \code{file_list}. Used to assign each sample to a group. The two groups represent the levels of the binary variable of interest.
}
  \item{seq_col}{
    Specifies the column of each sample's data frame containing the TCR/BCR sequences. Accepts either the column name or column index.
}
  \item{assoc_seqs}{
    A character vector containing the TCR/BCR sequences associated with the binary variable of interest.
}
  \item{nbd_radius}{
    The maximum distance (based on \code{dist_type}) between an associated sequence and other TCR/BCR sequences belonging to its neighborhood. Lower values require sequences to be more similar to an associated sequence in order to belong to its neighborhood.
}
  \item{dist_type}{
    Specifies the function used to quantify the similarity between sequences. The similarity between two sequences determines their pairwise distance, with greater similarity corresponding to shorter distance. Valid options are \code{"hamming"} (the default), which uses \code{\link{hamDistBounded}}, and \code{"levenshtein"}, which uses \code{\link{levDistBounded}}.
}
  \item{min_seq_length}{
    Clones with TCR/BCR sequences below this length will be removed. Passed to \code{\link{filterInputData}} when loading each sample.
}
  \item{drop_matches}{
    Passed to \code{\link{filterInputData}}. Accepts a regular expression (see \code{\link[base]{regex}}). Checks TCR/BCR sequences for a pattern match using \code{\link{grep}}. Those returning a match are dropped. By default, sequences containing either of the characters \code{*} or \code{_} are dropped.
}
  \item{subset_cols}{
    Controls which columns of the AIRR-Seq data from each sample are included in the output. Accepts a vector of column names or a vector of column indices. The default \code{NULL} includes all columns. Passed to \code{\link{filterInputData}}.
}
  \item{output_dir}{
  A file path to a directory for saving the output. A valid output directory is required, since no output is returned in R. The specified directory will be created if it does not already exist.
}
  \item{output_type}{
  A character string specifying the file format to use for saving the output. Valid options include \code{"csv"}, \code{"tsv"}, \code{"rds"} and \code{"rda"}.
}
  \item{verbose}{
  Logical. If \code{TRUE}, additional console output will be printed reporting the number of clones in each sample belonging to the neighborhood for each associated sequence, as well as the total number of clones in each neighborhood across all samples.
}
}
\details{
   For each associated sequence, its neighborhood is defined to include all clones with TCR/BCR sequences that are sufficiently similar to the associated sequence. The arguments \code{dist_type} and \code{nbd_radius} control how the similarity is measured and the degree of similarity required for neighborhood membership.

   For each associated sequence, a data frame is saved to an individual file. The data frame contains one row for each clone in the associated sequence's neighborhood (from all samples). It includes variables for sample ID, group ID and (if provided) subject ID, as well as variables from the AIRR-Seq data. Each neighborhood's file is named as its associated sequence, followed by the file extension based on \code{output_type}.

   The files saved by this function are intended for use with \code{\link{buildAssociatedClusterNetwork}}. See the \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters} package vignette for more details. The vignette can be accessed offline using \code{vignette("associated_clusters")}.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{findAssociatedSeqs}}
\code{\link{buildAssociatedClusterNetwork}}
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


## Step 2: Find Associated Clones ##

# output directory for current step
dir_nbds <- file.path(data_dir, "assoc_seq_nbds")

# Identify clones in a neighborhood around each associated sequence
findAssociatedClones(
  file_list = input_files, input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL, drop_matches = NULL,
  output_dir = dir_nbds)
}
}

