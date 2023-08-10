\name{buildPublicClusterNetwork}
\alias{buildPublicClusterNetwork}
\title{
    Build Global Network of Public TCR/BCR Clusters
}
\description{
    Part of the workflow \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}. Intended for use following \code{\link{findPublicClusters}}.

    Given node-level metadata for each sample's filtered clusters, combines the data into a global network and performs network analysis and cluster analysis.
}
\usage{
buildPublicClusterNetwork(

  ## Input ##
  file_list =
    list.files(file.path(getwd(), "public_clusters", "node_meta_data")),
  input_type = "rds",
  data_symbols = "ndat",
  header = TRUE, sep = "",
  seq_col,

  ## Network Settings ##
  drop_isolated_nodes = FALSE,
  node_stats = NULL, stats_to_include = NULL, cluster_stats = NULL,

  ## Visualization ##
  color_nodes_by = "SampleID",
  color_scheme = "turbo",
  plot_title = "Global Network of Public Clusters",

  ## Output ##
  output_dir =
    file.path(getwd(), "public_clusters"),
  output_name = "PublicClusterNetwork",

  ...

)
}
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the node-level metadata for one sample's filtered clusters. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{input_type}{
  A character string specifying the file format of the input files. Options are \code{"csv"}, \code{"rds"} and \code{"rda"}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of the data frame within each Rdata file. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
  Used for certain values of \code{input_type}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{seq_col}{
    Specifies the column in the node-level metadata that contains the TCR/BCR sequences. Accepts either the column name or column index.
}
  \item{drop_isolated_nodes}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{node_stats}{
    Deprecated. Does nothing since all network properties are automatically computed.
}
  \item{stats_to_include}{
    Deprecated. Does nothing since all network properties are automatically computed.
}
  \item{cluster_stats}{
    Deprecated. Does nothing since all network properties are automatically computed.
}
  \item{color_nodes_by}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{color_scheme}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{plot_title}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the global network.
}
  \item{output_dir}{
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
    The node-level metadata for the filtered clusters from all samples is combined in order to construct a single global network. Cluster analysis is used to partition the global network into clusters, which are considered as the public TCR/BCR clusters. Network properties for the nodes and clusters are computed and returned as metadata. A plot of the global network graph is produced, with the nodes colored according to sample ID.

    The computed node-level network properties are renamed to reflect their correspondence to the global network. This is done to distinguish them from the network properties that correspond to the sample-level networks. The names are \code{ClusterIDPublic}, \code{PublicNetworkDegree}, \code{PublicTransitivity}, \code{PublicCloseness}, \code{PublicCentralityByCloseness}, \code{PublicEigenCentrality}, \code{PublicCentralityByEigen}, \code{PublicBetweenness}, \code{PublicCentralityByBetweenness}, \code{PublicAuthorityScore}, \code{PublicCoreness}, \code{PublicPageRank}.

    See the \href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters} package vignette for more details. The vignette can be accessed offline using \code{vignette("public_clusters")}.
}
\value{
    A list of network objects as returned by \code{\link{buildRepSeqNetwork}}.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{findPublicClusters}}

\code{\link{buildPublicClusterNetworkByRepresentative}}
}
\examples{

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


## 2. Build Global Network of Public Clusters

# Directory of filtered node data
dir_filtered_samples_node <-
  file.path(dir_filtered_samples, "node_meta_data")

# Vector of file paths for filtered node data
files_filtered_samples_node <-
  list.files(dir_filtered_samples_node, full.names = TRUE)

dir_out <- file.path(data_dir, "public_clusters")

# Collect clones from all public clusters and perform network analysis
public_clusters <- buildPublicClusterNetwork(
  file_list = files_filtered_samples_node,
  seq_col = "CloneSeq", count_col = "CloneCount",
  size_nodes_by = 1,
  output_dir = dir_out)


# clean up temp directory
unlink(list.dirs(data_dir, full.names = TRUE),
       recursive = TRUE)


}
