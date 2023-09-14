\name{buildPublicClusterNetwork}
\alias{buildPublicClusterNetwork}
\title{
    Build Global Network of Public TCR/BCR Clusters
}
\description{
Part of the workflow
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}.
Intended for use following \code{\link[=findPublicClusters]{findPublicClusters()}}.

Given node-level metadata for each sample's filtered clusters, combines the
data into a global network and performs network analysis and cluster analysis.
}
\usage{
buildPublicClusterNetwork(

  ## Input ##
  file_list,
  input_type = "rds",
  data_symbols = "ndat",
  header = TRUE, sep,
  read.args = list(row.names = 1),
  seq_col,

  ## Network Settings ##
  drop_isolated_nodes = FALSE,
  node_stats = deprecated(),
  stats_to_include = deprecated(),
  cluster_stats = deprecated(),

  ## Visualization ##
  color_nodes_by = "SampleID",
  color_scheme = "turbo",
  plot_title = "Global Network of Public Clusters",

  ## Output ##
  output_dir = NULL,
  output_name = "PublicClusterNetwork",
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
  \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{input_type}{
  A character string specifying the file format of the input files. Options are
  \code{"csv"}, \code{"rds"} and \code{"rda"}. Passed to
  \code{\link[=loadDataFromFileList]{loadDataFromFileList()}}.
}
  \item{data_symbols}{
  Used when \code{input_type = "rda"}. Specifies the name of the data frame within
  each Rdata file. Passed to
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
  \item{seq_col}{
    Specifies the column in the node-level metadata that contains the TCR/BCR
    sequences.
    Accepts a character string containing the column name
      or a numeric scalar containing the column index.
}
  \item{drop_isolated_nodes}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{node_stats}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
   All network properties are automatically computed.
}
  \item{stats_to_include}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
  All network properties are automatically computed.
}
  \item{cluster_stats}{
  \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
  All network properties are automatically computed.
}
  \item{color_nodes_by}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
  The node-level network properties for the global network (see details) are
  included among the valid options.
}
  \item{color_scheme}{
  Passed to \code{\link[=addPlots]{addPlots()}}
  when constructing the global network.
}
  \item{plot_title}{
  Passed to \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}}
  when constructing the global network.
}
  \item{output_dir}{
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
  (including arguments to \code{\link[=addPlots]{addPlots()}})
  when constructing the global network.
  Does not include
  \code{node_stats}, \code{stats_to_include}, \code{cluster_stats}
  or \code{cluster_id_name}.
}
}
\details{
The node-level metadata for the filtered clusters from all samples is combined
and the global network is constructed by calling
\code{\link[=buildRepSeqNetwork]{buildNet()}} with
\code{node_stats = TRUE}, \code{stats_to_include = "all"},
\code{cluster_stats = TRUE} and \code{cluster_id_name = "ClusterIDPublic"}.

The computed node-level network properties are renamed to reflect their
correspondence to the global network. This is done to distinguish them from
the network properties that correspond to the sample-level networks. The names
are:

\itemize{
\item \code{ClusterIDPublic}
\item \code{PublicNetworkDegree}
\item \code{PublicTransitivity}
\item \code{PublicCloseness}
\item \code{PublicCentralityByCloseness}
\item \code{PublicEigenCentrality}
\item \code{PublicCentralityByEigen}
\item \code{PublicBetweenness}
\item \code{PublicCentralityByBetweenness}
\item \code{PublicAuthorityScore}
\item \code{PublicCoreness}
\item \code{PublicPageRank}
}

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters}
article on the package website.
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

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html}{Searching for Public TCR/BCR Clusters article on package website}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link[=findPublicClusters]{findPublicClusters()}}

\code{\link[=buildPublicClusterNetworkByRepresentative]{buildPublicClusterNetworkByRepresentative()}}
}
\examples{
set.seed(42)

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
    plot_title = NULL,
    plot_subtitle = NULL,
    print_plots = TRUE
  )


\dontshow{
# Clean up temporary files
file.remove(
  file.path(tempdir(),
            paste0("Sample", 1:samples, ".rds")
  )
)
unlink(
  file.path(tempdir(), c("node_meta_data", "cluster_meta_data")),
  recursive = TRUE
)
}
}
