\name{buildAssociatedClusterNetwork}
\alias{buildAssociatedClusterNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Combine Associated Sequence Neighborhoods and Build Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    After running \code{findAssociatedClones}, this function is used to combine
    some or all of the resulting neighborhoods into a single network in which
    clustering and network analysis are performed.
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{input_type}{
%%     ~~Describe \code{input_type} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{data_symbols}{
%%     ~~Describe \code{data_symbols} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{node_stats}{
%%     ~~Describe \code{node_stats} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{cluster_stats}{
%%     ~~Describe \code{cluster_stats} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{output_name}{
%%     ~~Describe \code{output_name} here~~
  Passed to \code{buildRepSeqNetwork}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
  Other arguments to \code{buildRepSeqNetwork}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    Essentially a wrapper to \code{loadDataFromFileList} and
    \code{buildRepSeqNetwork}, with default argument values tailored to the
    associated cluster workflow.
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
\code{\link{findAssociatedSeqs}}
\code{\link{findAssociatedClones}}
}
\examples{
\dontrun{
## Generate some toy data for demonstration

# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_input_samples <- file.path(data_dir, "input_samples")
dir.create(dir_input_samples, showWarnings = FALSE)

samples <- 30
affixes <- c("AAAA", "AASA", "AACA", "AAQA", "AAQ", "AAA", "AASAA", "AAAAA")
affix_probs_g0 <- rep(1 / length(affixes),
                      times = length(affixes) * samples / 2)
affix_probs_g1 <- rep(c(1, 5, 1, 1, 1, 1, 5, 1), times = samples / 2)
affix_probs <- matrix(c(affix_probs_g0, affix_probs_g1),
                      nrow = samples, byrow = TRUE)
new_probs_g0 <- rep(c(1/2, 1/6, 1/6, 1/6), times = samples / 2)
new_probs_g1 <- rep(c(1/3, 1/6, 1/6, 1/3), times = samples / 2)
new_probs <- matrix(c(new_probs_g0, new_probs_g1),
                    nrow = samples, byrow = TRUE)
simulateToyData(
  samples = samples,
  sample_size = 30,
  prefix_length = 1,
  prefix_chars = c("A", "C"),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = affixes,
  affix_probs = affix_probs,
  num_edits = 4,
  edit_pos_probs = function(seq_length) {
    dnorm(seq(-4, 4, length.out = seq_length))
  },
  edit_ops = c("insertion", "deletion", "transmutation"),
  edit_probs = c(5, 1, 4),
  new_chars = c("A", "S", "C", "Q"),
  new_probs = new_probs,
  output_dir = dir_input_samples,
  no_return = TRUE
)


## 1. Find Associated Sequences

# input files for step 1 (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))
head(input_files)

# group labels for the samples
group_labels <- c(rep("reference", samples / 2), rep("comparison", samples / 2))

# search across samples for associated sequences using Fisher's exact test
associated_seqs <- findAssociatedSeqs(
  file_list = input_files, input_type = "rds",
  group_ids = group_labels, groups = c("reference", "comparison"),
  min_seq_length = NULL, drop_matches = NULL,
  seq_col = "CloneSeq", outfile = NULL)
head(associated_seqs)


## 2. Find Associated Clones

# output directory for current step
dir_nbds <- file.path(data_dir, "assoc_seq_nbds")

# Identify neighborhood around each associated sequence
findAssociatedClones(
  file_list = input_files, input_type = "rds", group_ids = group_labels,
  seq_col = "CloneSeq", dist_type = "levenshtein",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL, drop_matches = NULL,
  output_dir = dir_nbds)


## 3. Build Associated Cluster Network

# Files created during previous step
nbd_files <- list.files(dir_nbds, full.names = TRUE)
nbd_files

# Combine neighborhoods and perform network analysis
all_clusters <- buildAssociatedClusterNetwork(
  file_list = nbd_files,
  seq_col = "CloneSeq", dist_type = "levenshtein", size_nodes_by = 1.5,
  output_dir = file.path(data_dir, "assoc_clusters"))

# focus on a particular cluster
buildRepSeqNetwork(
  data = all_clusters$node_data[all_clusters$node_data$cluster_id == 3, ],
  seq_col = "CloneSeq", color_nodes_by = "CloneSeq", size_nodes_by = 2,
  output_dir = NULL, output_name = "Cluster 3")
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
