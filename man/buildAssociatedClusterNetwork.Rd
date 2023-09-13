\name{buildAssociatedClusterNetwork}
\alias{buildAssociatedClusterNetwork}
\title{
    Build Global Network of Associated TCR/BCR Clusters
}
\description{
Part of the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}.
Intended for use following
\code{\link[=findAssociatedClones]{findAssociatedClones()}}.

Given data containing a neighborhood of similar clones around each associated sequence,
combines the data into a global network and performs network analysis and cluster
analysis.
}
\usage{
buildAssociatedClusterNetwork(
  file_list,
  input_type = "rds",
  data_symbols = "data", header = TRUE, sep,
  read.args = list(row.names = 1),
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
  verbose = FALSE,
  ...
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
  A character string specifying the file format of the neighborhood data files.
  Options are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"},
  \code{"rds"} and \code{"rda"}.
  Passed to \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each neighborhood's
  data frame within its respective Rdata file.  Passed to
  \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{header}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument is used to specify the value of the \code{header}
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
    this argument is used to specify values of optional
    arguments to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
    Accepts a named list of argument values.
    Values of \code{header} and \code{sep}
    in this list take precedence over values specified via the \code{header}
    and \code{sep} arguments.
}
  \item{seq_col}{
    Specifies the column of each neighborhood's data frame containing the
    TCR/BCR sequences.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
  \item{min_seq_length}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{drop_matches}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{drop_isolated_nodes}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{node_stats}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{stats_to_include}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{cluster_stats}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{color_nodes_by}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{output_name}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
  \item{\dots}{
  Other arguments to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
}
\details{
Each associated sequence's neighborhood contains clones (from all samples)
with TCR/BCR sequences similar to the associated sequence. The neighborhoods
are assumed to have been previously identified using
\code{\link[=findAssociatedClones]{findAssociatedClones()}}.

The neighborhood data for all associated sequences are used to construct a
single global network. Cluster analysis is used to partition the global network
into clusters, which are considered as the associated TCR/BCR clusters. Network
properties for the nodes and clusters are computed and returned as metadata.
A plot of the global network graph is produced, with the nodes colored according
to the binary variable of interest.

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html}{Searching for Associated TCR/BCR Clusters}
article on the package website for more details.
}
\value{
    A list of network objects as returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}.
    The list is returned invisibly.
    If the input data contains a combined total of fewer than two rows, or if the
    global network contains no nodes, then the function returns \code{NULL},
    invisibly, with a warning.
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
\code{\link[=findAssociatedClones]{findAssociatedClones()}}
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