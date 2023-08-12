\name{buildPublicClusterNetworkByRepresentative}
\alias{buildPublicClusterNetworkByRepresentative}
\title{
    Build Global Network of Public TCR/BCR Clusters Using Representative Clones
}
\description{
Alternative step in the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}.
Intended for use following \code{\link{findPublicClusters}} in cases where
\code{\link{buildPublicClusterNetwork}} cannot be practically used due to the
size of the full global network.

Given cluster-level metadata for each sample's filtered clusters, selects a
representative TCR/BCR from each cluster, combines the representatives into a
global network and performs network analysis and cluster analysis.
}
\usage{
buildPublicClusterNetworkByRepresentative(

  ## Input ##
  file_list,
  input_type = "rds",
  data_symbols = "cdat",
  header = TRUE, sep = "",
  seq_col = "seq_w_max_count",
  count_col = "agg_count",

  ## Network Settings ##
  dist_type = "hamming",
  dist_cutoff = 1,
  cluster_fun = cluster_fast_greedy,

  ## Visualization ##
  plots = TRUE,
  print_plots = TRUE,
  plot_title = "auto",
  plot_subtitle = "auto",
  color_nodes_by = "SampleID",
  color_scheme = "turbo",
  ...,

  ## Output ##
  output_dir =
    file.path(getwd(), "public_clusters"),
  output_type = "rda",
  output_name = "PubClustByRepresentative",
  pdf_width = 12,
  pdf_height = 10

)
}
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the cluster-level metadata for
  one sample's filtered clusters. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{input_type}{
  A character string specifying the file format of the input files. Options are
  \code{"csv"}, \code{"rds"} and \code{"rda"}. Passed to
  \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of the data frame
  within each Rdata file. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
  Used for certain values of \code{input_type}. Passed to
  \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
  Used for certain values of \code{input_type}. Passed to
  \code{\link{loadDataFromFileList}}.
}
  \item{seq_col}{
    Specifies the column in the cluster-level metadata that contains the
    representative TCR/BCR sequence for each cluster. Accepts either the column
    name or column index. By default, uses the sequence with the maximum clone
    count in each cluster.
}
  \item{count_col}{
    Specifies the column in the cluster-level metadata that contains the aggregate
    clone count for each cluster. Accepts either the column name or column index.
}
  \item{dist_type}{
    Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{dist_cutoff}{
    Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{cluster_fun}{
    Passed to \code{\link{buildRepSeqNetwork}} when performing cluster analysis on
    the global network.
}
  \item{plots}{
    Logical. Should plots of the global network graph be produced?
}
  \item{print_plots}{
    Logical. If plots of the global network graph are produced, should they be
    printed to the R plotting window?
}
  \item{plot_title}{
    Passed to \code{\link{generateNetworkGraphPlots}} when producing plots of the
    global network graph.
}
  \item{plot_subtitle}{
    Passed to \code{\link{generateNetworkGraphPlots}} when producing plots of the
    global network graph.
}
  \item{color_nodes_by}{
    Passed to \code{\link{generateNetworkGraphPlots}} when producing plots of the
    global network graph.
}
  \item{color_scheme}{
    Passed to \code{\link{generateNetworkGraphPlots}} when producing plots of the
    global network graph.
}
  \item{\dots}{
    Other arguments to \code{\link{generateNetworkGraphPlots}} when producing plots
    of the global network graph.
}
  \item{output_dir}{
    Passed to \code{saveNetwork} after constructing the global network.
}
  \item{output_type}{
    Passed to \code{saveNetwork} after constructing the global network.
}
  \item{output_name}{
    Passed to \code{saveNetwork} after constructing the global network.
}
  \item{pdf_width}{
    Passed to \code{saveNetwork} after constructing the global network.
    Only applicable if \code{plots = TRUE}.
}
  \item{pdf_height}{
    Passed to \code{saveNetwork} after constructing the global network.
    Only applicable if \code{plots = TRUE}.
}
}
\details{
From each filtered cluster in each sample's network, a representative TCR/BCR
is selected. By default, this is the sequence with the greatest clone count in
each cluster. The representatives from all clusters and all samples are then
used to construct a single global network. Cluster analysis is used to partition
this global network into clusters. Network properties for the nodes and clusters
are computed and returned as metadata. A plot of the global network graph is
produced, with the nodes colored according to sample ID.

Within this network, clusters containing nodes from multiple samples can be
considered as the skeletons of the complete public clusters. The filtered cluster
data for each sample can then be subset to keep the sample-level clusters whose
representative TCR/BCRs belong to the skeletons of the public clusters. After
subsetting in this manner, \code{\link{buildPublicClusterNetwork}} can be used
to construct the global network of complete public clusters.

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}
package vignette for more details. The vignette can be accessed offline using
\code{vignette("public_clusters")}.
}
\value{
    If the input data contains a combined total of fewer than two rows, or if the
    global network contains no nodes, then the function returns \code{NULL},
    invisibly, with a warning. Otherwise, invisibly returns
    a list of network objects as returned by \code{\link{buildRepSeqNetwork}}.
    The data frame \code{cluster_data} includes the following variables that
    represent properties of the clusters in the global network of representative
    TCR/BCRs:
\item{TotalSampleLevelNodes}{
    For each representative TCR/BCR in the global cluster, we record the number
    of nodes in the sample-level cluster for which it is the representative TCR/BCR.
    We then sum these node counts across all the representative TCR/BCRs in the
    global cluster.}
\item{TotalCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the
    aggregate clone count from all nodes in the sample-level cluster for which
    it is the representative TCR/BCR. We then sum these aggregate clone counts
    across all the representative TCR/BCRs in the global cluster.}
\item{MeanOfMeanSeqLength}{
    For each representative TCR/BCR in the global cluster, we record the mean
    sequence length over all clones (nodes) in the sample-level cluster for which
    it is the representative TCR/BCR. We then average these mean sequence lengths
    over all the representative TCR/BCRs in the global cluster.}
\item{MeanDegreeInPublicNet}{
    For each representative TCR/BCR in the global cluster, we record the mean
    network degree over all nodes in the sample-level cluster for which it is the
    representative TCR/BCR. We then average these mean degree values over all the
    representative TCR/BCRs in the global cluster.}
\item{MaxDegreeInPublicNet}{
    For each representative TCR/BCR in the global cluster, we record the maximum
    network degree across all nodes in the sample-level cluster for which it is
    the representative TCR/BCR. We then take the maximum of these maximum degree
    values over all the representative TCR/BCRs in the global cluster.}
\item{SeqWithMaxDegree}{
    For each representative TCR/BCR in the global cluster, we record the maximum
    network degree across all nodes in the sample-level cluster for which it is
    the representative TCR/BCR. We then identify the representative TCR/BCR with
    the maximum value of these maximum degrees over all the representative TCR/BCRs
    in the global cluster. The TCR/BCR sequence of the identified representative
    TCR/BCR is recorded in this variable.}
\item{MaxCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the maximum
    clone count across all clones (nodes) in the sample-level cluster for which
    it is the representative TCR/BCR. We then take the maximum of these maximum
    clone counts over all the representative TCR/BCRs in the global cluster.}
\item{SampleWithMaxCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the maximum
    clone count across all clones (nodes) in the sample-level cluster for which
    it is the representative TCR/BCR. We then identify the representative TCR/BCR
    with the maximum value of these maximum clone counts over all the representative
    TCR/BCRs in the global cluster. The sample to which the identified representative
    TCR/BCR belongs is recorded in this variable.}
\item{SeqWithMaxCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the maximum
    clone count across all clones (nodes) in the sample-level cluster for which it
    is the representative TCR/BCR. We then identify the representative TCR/BCR with
    the maximum value of these maximum clone counts over all the representative
    TCR/BCRs in the global cluster. The TCR/BCR sequence of the identified
    representative TCR/BCR is recorded in this variable.}
\item{MaxAggCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the
    aggregate clone count across all clones (nodes) in the sample-level cluster
    for which it is the representative TCR/BCR. We then take the maximum of these
    aggregate clone counts over all the representative TCR/BCRs in the global
    cluster.}
\item{SampleWithMaxAggCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the
    aggregate clone count across all clones (nodes) in the sample-level cluster
    for which it is the representative TCR/BCR. We then identify the representative
    TCR/BCR with the maximum value of these aggregate clone counts over all the
    representative TCR/BCRs in the global cluster. The sample to which the
    identified representative TCR/BCR belongs is recorded in this variable.}
\item{SeqWithMaxAggCloneCount}{
    For each representative TCR/BCR in the global cluster, we record the
    aggregate clone count across all clones (nodes) in the sample-level cluster
    for which it is the representative TCR/BCR. We then identify the representative
    TCR/BCR with the maximum value of these aggregate clone counts over all the
    representative TCR/BCRs in the global cluster. The TCR/BCR sequence of the
    identified representative TCR/BCR is recorded in this variable.}
\item{DiameterLength}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{Assortativity}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{GlobalTransitivity}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{EdgeDensity}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{DegreeCentralityIndex}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{ClosenessCentralityIndex}{
  See \code{\link{getClusterStats}}.
  Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{EigenCentralityIndex}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}
\item{EigenCentralityEigenvalue}{
    See \code{\link{getClusterStats}}.
    Based on edge connections between representative TCR/BCRs in the global cluster.
}


}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{findPublicClusters}}
\code{\link{buildPublicClusterNetwork}}
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

# relative generation probabilities
pgen <- cbind(
  stats::toeplitz(0.6^(0:(sample_size - 1))),
  matrix(1, nrow = samples, ncol = length(base_seqs) - samples))

# Simulate the data
foo <- simulateToyData(
  samples = samples, sample_size = sample_size,
  prefix_length = 1, prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = base_seqs, affix_probs = pgen, num_edits = 0,
  output_dir = dir_input_samples, no_return = FALSE)

# View first few rows of data for sample 1
head(readRDS(file.path(dir_input_samples, "Sample1.rds")))


## 1. Find Public Clusters in Each Sample

# input files for step 1 (one per sample)
input_files <-
  file.path(dir_input_samples,
            paste0("Sample", 1:samples, ".rds"))
head(input_files)

# directory to save filtered sample data
dir_filtered_samples <- file.path(data_dir, "filtered_samples")

# Filter Clusters Within Each Sample
findPublicClusters(
  file_list = input_files, input_type = "rds",
  seq_col = "CloneSeq", count_col = "CloneCount",
  min_seq_length = NULL, drop_matches = NULL,
  top_n_clusters = 3, min_node_count = 5, min_clone_count = 15000,
  output_dir = dir_filtered_samples)


## 2. Build Public Cluster Network by Representative TCR/BCRs

# Cluster-level meta data for each sample's public clusters
dir_filtered_samples_cluster <-
  file.path(dir_filtered_samples, "cluster_meta_data")
files_filtered_samples_cluster <-
  list.files(dir_filtered_samples_cluster, full.names = TRUE)
head(files_filtered_samples_cluster)

dir_out <- file.path(data_dir, "public_clusters")

buildPublicClusterNetworkByRepresentative(
  file_list = files_filtered_samples_cluster,
  output_dir = dir_out,
  size_nodes_by = 1)


# clean up temp directory
unlink(list.dirs(data_dir, full.names = TRUE),
       recursive = TRUE)

}
}