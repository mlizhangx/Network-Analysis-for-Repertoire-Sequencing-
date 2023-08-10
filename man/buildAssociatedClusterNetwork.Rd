\name{buildAssociatedClusterNetwork}
\alias{buildAssociatedClusterNetwork}
\title{
    Build Global Network of Associated TCR/BCR Clusters
}
\description{
    Part of the workflow \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}. Intended for use following \code{\link{findAssociatedClones}}.

    Given data containing a neighborhood of similar clones around each associated sequence, combines the data into a global network and performs network analysis and cluster analysis.
}
\usage{
buildAssociatedClusterNetwork(
    file_list,
    input_type = "csv",
    data_symbols = "data",
    header = TRUE, sep = ",",
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE,
    stats_to_include =
      chooseNodeStats(cluster_id = TRUE),
    cluster_stats = TRUE,
    color_nodes_by = "GroupID",
    output_name = "AssociatedClusterNetwork",
    ...
)
}
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the neighborhood data for an associated sequence. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{input_type}{
  A character string specifying the file format of the neighborhood data files. Options are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"}, \code{"rds"} and \code{"rda"}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each neighborhood's data frame within its respective Rdata file.  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{seq_col}{
    Specifies the column of each neighborhood's data frame containing the TCR/BCR sequences. Accepts either the column name or column index.
}
  \item{min_seq_length}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{drop_matches}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{drop_isolated_nodes}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{node_stats}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{stats_to_include}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{cluster_stats}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{color_nodes_by}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{output_name}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{\dots}{
  Other arguments to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
}
\details{
    Each associated sequence's neighborhood contains clones (from all samples) with TCR/BCR sequences similar to the associated sequence. The neighborhoods are assumed to have been previously identified using \code{\link{findAssociatedClones}}.

    The neighborhood data for all associated sequences are used to construct a single global network. Cluster analysis is used to partition the global network into clusters, which are considered as the associated TCR/BCR clusters. Network properties for the nodes and clusters are computed and returned as metadata. A plot of the global network graph is produced, with the nodes colored according to the binary variable of interest.

    See the \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters} package vignette for more details. The vignette can be accessed offline using \code{vignette("associated_clusters")}.
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
\code{\link{findAssociatedSeqs}}
\code{\link{findAssociatedClones}}
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


## Step 3: Global Network of Associated Clusters ##

# Files created during previous step
nbd_files <- list.files(dir_nbds, full.names = TRUE)
nbd_files

# Combine neighborhoods and perform network analysis
all_clusters <- buildAssociatedClusterNetwork(
  file_list = nbd_files,
  seq_col = "CloneSeq",
  size_nodes_by = 1.5,
  output_dir = file.path(data_dir, "assoc_clusters"))


# clean up temp directory
unlink(list.dirs(data_dir, full.names = TRUE),
       recursive = TRUE)


}