\name{findPublicClusters}
\alias{findPublicClusters}
\title{
    Find Public Clusters Among RepSeq Samples
}
\description{
Part of the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}.

Given multiple samples of bulk Adaptive Immune Receptor Repertoire Sequencing
(AIRR-Seq) data, construct the repertoire network for each sample. Within each
sample's network, perform cluster analysis and filter the clusters based on node
count and aggregate clone count.
}
\usage{
findPublicClusters(

  ## Input ##
  file_list,
  input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids =
    paste0("Sample", 1:length(file_list)),
  seq_col,
  count_col = NULL,

  ## Search Criteria ##
  min_seq_length = 3,
  drop_matches = "[*|_]",
  top_n_clusters = 20,
  min_node_count = 10,
  min_clone_count = 100,

  ## Optional Visualization ##
  plots = FALSE,
  print_plots = FALSE,
  plot_title = "auto",
  color_nodes_by = "cluster_id",

  ## Output ##
  output_dir =
    file.path(getwd(), "public_clusters"),
  output_type = "rds",

  ## Optional Output ##
  output_dir_unfiltered = NULL,
  output_type_unfiltered = "rds",

  ...

)
}
\arguments{
  \item{file_list}{
  A vector of file paths where each file contains the data for a sample, with
  observations indexed by row. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{input_type}{
  A character string specifying the file format of the sample data files. Options
  are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"}, \code{"rds"} and
  \code{"rda"}. Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of each sample's data
  frame within its respective Rdata file.
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
  Used for certain values of \code{input_type}.
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
  Used for certain values of \code{input_type}.
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sample_ids}{
  A character or numeric vector of sample IDs, whose length matches that of
  \code{file_list}. The values should be valid for use as filenames and should
  avoid using the forward slash or backslash characters (\code{/} or \code{\\}).
}
  \item{seq_col}{
    Specifies the column of each sample's data frame containing the TCR/BCR
    sequences. Accepts either the column name or column index.
}
  \item{count_col}{
    Specifies the column of each sample's data frame containing the clone count
    (measure of clonal abundance). Accepts either the column name or column index.
    If \code{NULL}, the clusters in each sample's network will be selected solely
    based upon node count.
}
  \item{min_seq_length}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample.
}
  \item{drop_matches}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample. Accepts a character string containing a regular expression
  (see \code{\link[base]{regex}}).
  Checks TCR/BCR sequences for a pattern match. Those returning a match are dropped.
  By default, sequences containing any of the characters \code{*}, \code{|}
  or \code{_} are dropped.
}
  \item{top_n_clusters}{
    The number of clusters from each sample to be automatically be included among
    the filtered clusters, based on greatest node count.
}
  \item{min_node_count}{
    Clusters with at least this many nodes will be included among the filtered
    clusters.
}
  \item{min_clone_count}{
    Clusters with an aggregate clone count of at least this value will be included
    among the filtered clusters. A value of \code{NULL} ignores this criterion
    and does not select additional clusters based on clone count.
}
  \item{plots}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample.
}
  \item{print_plots}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample.
}
  \item{plot_title}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample.
}
  \item{color_nodes_by}{
  Passed to \code{\link{buildRepSeqNetwork}} when constructing the network for
  each sample.
}
  \item{output_dir}{
    The file path of the directory for saving the output. The directory will be
    created if it does not already exist.
}
  \item{output_type}{
  A character string specifying the file format to use for saving the output.
  Valid options include \code{"csv"}, \code{"rds"} and \code{"rda"}.
}
  \item{output_dir_unfiltered}{
    An optional directory for saving the unfiltered network data for each sample.
    By default, only the filtered results are saved.
}
  \item{output_type_unfiltered}{
    A character string specifying the file format to use for saving the unfiltered
    network data for each sample. Only applicable if \code{output_dir_unfiltered}
    is non-null. Passed to \code{\link{buildRepSeqNetwork}} when constructing the
    network for each sample.
}
  \item{\dots}{
    Other arguments to \code{buildRepSeqNetwork} when constructing the network
    for each sample, not including \code{node_stats}, \code{stats_to_include} or
    \code{cluster_stats} (see details).
}
}
\details{
All node-level and cluster-level network properties are automatically computed for
each sample's network. The node-level properties are renamed to reflect their
correspondence to the sample-level network. Specifically, the properties are named
\code{ClusterIDInSample}, \code{SampleLevelNetworkDegree},
\code{SampleLevelTransitivity}, \code{SampleLevelCloseness},
\code{SampleLevelCentralityByCloseness}, \code{SampleLevelCentralityByEigen},
\code{SampleLevelEigenCentrality}, \code{SampleLevelBetweenness},
\code{SampleLevelCentralityByBetweenness}, \code{SampleLevelAuthorityScore},
\code{SampleLevelCoreness}, \code{SampleLevelPageRank}. A variable \code{SampleID}
is added to both the node-level and cluster-level meta data for each sample.

After the clusters in each sample are filtered, the node-level and cluster-level
metadata are saved in the respective subdirectories \code{node_meta_data} and
\code{cluster_meta_data} of the output directory specified by \code{output_dir};
each is saved using the sample ID (the corresponding element of
\code{sample_id_list}) as a filename prefix, followed by the appropriate file
extension based on \code{output_type}.

The unfiltered network results for each sample can also be saved by supplying a
directory to \code{output_dir_unfiltered}, if these results are desired for
downstream analysis.

The files containing the node-level metadata for the filtered clusters can be
supplied to \code{\link{buildPublicClusterNetwork}} in order to construct a global
network of public clusters. If the full global network is too large to practically
construct, the files containing the cluster-level meta data for the filtered
clusters can be supplied to \code{\link{buildPublicClusterNetworkByRepresentative}}
to build a global network using only a single representative sequence from each
cluster. This allows prominent public clusters to still be identified.

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}
package vignette for more details. The vignette can be accessed offline using
\code{vignette("public_clusters")}.

}
\value{
    Returns \code{TRUE}, invisibly.
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
\code{\link{buildPublicClusterNetwork}}

\code{\link{buildPublicClusterNetworkByRepresentative}}
}
\examples{

## Simulate 30 samples with a mix of public/private sequences ##
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
# Relative generation probabilities
pgen <- cbind(
  stats::toeplitz(0.6^(0:(sample_size - 1))),
  matrix(1, nrow = samples, ncol = length(base_seqs) - samples)
)
simulateToyData(
  samples = samples,
  sample_size = sample_size,
  prefix_length = 1,
  prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = base_seqs,
  affix_probs = pgen,
  num_edits = 0,
  output_dir = tempdir(),
  no_return = TRUE
)


## 1. Find Public Clusters in Each Sample
sample_files <-
  file.path(tempdir(),
            paste0("Sample", 1:samples, ".rds")
  )
findPublicClusters(
  file_list = sample_files,
  input_type = "rds",
  seq_col = "CloneSeq",
  count_col = "CloneCount",
  min_seq_length = NULL,
  drop_matches = NULL,
  top_n_clusters = 3,
  min_node_count = 5,
  min_clone_count = 15000,
  output_dir = tempdir()
)

## 2. Build Global Network of Public Clusters
public_clusters <-
  buildPublicClusterNetwork(
    file_list =
      list.files(
        file.path(tempdir(), "node_meta_data"),
        full.names = TRUE
      ),
    seq_col = "CloneSeq",
    count_col = "CloneCount",
    size_nodes_by = 1,
    output_dir = tempdir()
  )


# Clean up temporary files
file.remove(
  file.path(tempdir(),
            c(paste0("Sample", 1:samples, ".rds"),
              "PublicClusterNetwork_NodeMetadata.csv",
              "PublicClusterNetwork_ClusterMetadata.csv",
              "PublicClusterNetwork.pdf",
              "PublicClusterNetwork_EdgeList.txt",
              "PublicClusterNetwork_AdjacencyMatrix.mtx"
            )
  )
)
unlink(
  file.path(tempdir(), c("node_meta_data", "cluster_meta_data")),
  recursive = TRUE
)

}
