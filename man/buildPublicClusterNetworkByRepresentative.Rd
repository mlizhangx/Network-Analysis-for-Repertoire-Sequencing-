\name{buildPublicClusterNetworkByRepresentative}
\alias{buildPublicClusterNetworkByRepresentative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Build Public Cluster Network By Representative
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given cluster-level meta data for each sample's filtered public clusters,
    combine the data across samples and perform network analysis using a
    representative sequence from each cluster.
}
\usage{
buildPublicClusterNetworkByRepresentative(

  ## Input ##
  file_list,
  input_type = "rds", data_symbols = "cdat",
  header = TRUE, sep = "",
  seq_col = "seq_w_max_count",
  count_col = "agg_count",

  ## Network ##
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
  output_dir = file.path(getwd(), "public_clusters"),
  output_type = "rda",
  output_name = "PubClustByRepresentative",
  pdf_width = 12, pdf_height = 10

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} when loading the cluster-level meta data
  for each sample (one data frame per sample). Data frames should match the
  format of the cluster-level meta data returned by \code{buildRepSeqNetwork}.
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
    The column name or number of each data frame that contains the representative sequence to be used as the basis of similarity between rows during network analysis. This column must have the same name in each sample.
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{dist_cutoff}{
%%     ~~Describe \code{dist_cutoff} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{cluster_fun}{
%%     ~~Describe \code{net} here~~
    Passed to \code{addClusterMembership}. Controls the algorithm used for cluster identification.
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{print_plots}{
%%     ~~Describe \code{print_plots} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{plot_title}{
%%     ~~Describe \code{plot_title} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{plot_subtitle}{
%%     ~~Describe \code{plot_subtitle} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Passed to \code{generateNetworkGraphPlots}.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    Passed to \code{generateNetworkGraphPlots}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
    Other arguments to \code{generateNetworkGraphPlots}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    Passed to \code{saveNetwork}.
}
  \item{output_type}{
%%     ~~Describe \code{output_type} here~~
    Passed to \code{saveNetwork}.
}
  \item{output_name}{
%%     ~~Describe \code{output_name} here~~
    Passed to \code{saveNetwork}.
}
  \item{pdf_width}{
%%     ~~Describe \code{pdf_width} here~~
    Passed to \code{saveNetwork}.
}
  \item{pdf_height}{
%%     ~~Describe \code{pdf_height} here~~
    Passed to \code{saveNetwork}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    This function is intended for performing network analysis on the public clusters
    obtained using \code{findPublicClusters()}.

    By using the cluster-level meta data as the input to \code{buildRepSeqNetwork()},
    this function treats each public cluster as a single node and performs network
    analysis with similarity based on a representative sequence for each cluster
    (e.g., the sequence with the greatest clone count). \code{buildRepSeqNetwork} is
    called without any filtering and with \code{drop_isolated_nodes = FALSE}.

    All node-level properties are automatically computed for the network, and
    cluster-level properties for the network (in which clusters of nodes represent
    clusters of public-cluster-representatives) are computed based on the properties
    of the public cluster represented by each node. See the 'value' section for a list and      description of each property.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A list of network objects as per \code{buildRepSeqNetwork()}. The data frame
    \code{cluster_data} includes the following variables:
    \item{TotalSampleLevelNodes}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the sum of the sample-network-level nodes
        contained in all such public clusters for the current row.}
    \item{TotalCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the sum of the aggregate clone count
        across all such public clusters for the current row.}
    \item{MeanOfMeanSeqLength}{
        Each network node in the cluster for the current row represents
        a public cluster; this value is the mean value of the mean sequence length
        across all such public clusters for the current row.}
    \item{MeanDegreeInPublicNet}{
        The mean network degree of all the network nodes in the cluster for
        the current row.}
    \item{MaxDegreeInPublicNet}{
        The maximum network degree of all the network nodes in the cluster for
        the current row.}
    \item{SeqWithMaxDegree}{
        The representative sequence of the network node with maximum network degree
        of all the network nodes in the cluster for the current row. If more than
        one node attains the maximum network degree, the first sequence found is
        returned.}
    \item{MaxCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this is the maximum value of the
        maximum clone count property across all such public clusters for the
        current row.}
    \item{SampleWithMaxCloneCount}{
        The sample possessing the public cluster with the maximum value
        of the maximum clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sample ID found is returned.}
    \item{SeqWithMaxCloneCount}{
        The representative sequence of the public cluster with the maximum value
        of the maximum clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sequence found is returned.}
    \item{MaxAggCloneCount}{
        Each network node in the cluster for the current row represents
        a public cluster; this is the maximum value of the
        aggregate clone count property across all such public clusters for the
        current row.}
    \item{SampleWithMaxAggCloneCount}{
        The sample possessing the public cluster with the maximum value
        of the aggregate clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sample ID found is returned.}
    \item{SeqWithMaxAggCloneCount}{
        The representative sequence of the public cluster with the maximum value
        of the aggregate clone count property across all the public clusters
        represented by the network nodes in the cluster for the current row.
        If more than one public cluster attains the maximum value, the first
        sequence found is returned.}
    \item{DiameterLength}{
        See \code{?getClusterStats()}.}
    \item{Assortativity}{
        See \code{?getClusterStats()}.}
    \item{GlobalTransitivity}{
        See \code{?getClusterStats()}.}
    \item{EdgeDensity}{
        See \code{?getClusterStats()}.}
    \item{DegreeCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{ClosenessCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{EigenCentralityIndex}{
        See \code{?getClusterStats()}.}
    \item{EigenCentralityEigenvalue}{
        See \code{?getClusterStats()}.}
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


## 1. Find Public Clusters in Each Sample

# input files for step 1 (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))
head(input_files)

# Search across samples for public clusters
dir_filtered_samples <- file.path(data_dir, "filtered_samples")
findPublicClusters(
  file_list = input_files, input_type = "rds",
  sample_ids = paste0("Sample", 1:samples),
  seq_col = "CloneSeq", count_col = "CloneCount",
  min_seq_length = NULL, drop_matches = NULL,
  output_dir = dir_filtered_samples)


## 2. Build Public Cluster Network by Representative Sequence

# Cluster-level meta data for each sample's public clusters
dir_filtered_samples_cluster <- file.path(dir_filtered_samples, "cluster_meta_data")
files_filtered_samples_cluster <- list.files(dir_filtered_samples_cluster,
                                             full.names = TRUE)
head(files_filtered_samples_cluster)

dir_out <- file.path(data_dir, "public_clusters")

buildPublicClusterNetworkByRepresentative(
  file_list = files_filtered_samples_cluster,
  color_nodes_by = "ClusterIDPublic", output_dir = dir_out,
  size_nodes_by = 1)
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
