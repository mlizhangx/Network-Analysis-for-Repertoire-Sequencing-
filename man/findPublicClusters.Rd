\name{findPublicClusters}
\alias{findPublicClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Find Public Clusters Among RepSeq Samples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple samples of RepSeq data, perform network analysis on each
    sample individually and filter clusters by node count and clone count.
}
\usage{
findPublicClusters(

  ## Input ##
  file_list, input_type,
  data_symbols = NULL, header = TRUE, sep = "",
  sample_ids = 1:length(file_list),
  seq_col, count_col = NULL,
  min_seq_length = 3, drop_matches = "[*|_]",

  ## Network ##
  top_n_clusters = 20,
  min_node_count = 10,
  min_clone_count = 100,

  ## Visualization ##
  plots = FALSE, print_plots = FALSE,
  plot_title = "auto", color_nodes_by = "cluster_id",

  ## Output ##
  output_dir = file.path(getwd(), "public_clusters"),
  output_type = "rds",
  output_dir_unfiltered = NULL,
  output_type_unfiltered = "rds",
  ...

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList}(one data frame per sample).
}
  \item{input_type}{
%%     ~~Describe \code{input_type} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{data_symbols}{
%%     ~~Describe \code{data_symbols} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
  Passed to \code{\link{loadDataFromFileList}}.
}
  \item{sample_ids}{
%%     ~~Describe \code{sample_ids} here~~
  A character or numeric vector of sample IDs, whose length matches that of \code{file_list}. The values
  should be valid for use in filenames.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of each sample's data frame that contains the receptor sequences to be used as the basis of similarity between rows during network analysis. This column must have the same name in each sample.
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
    The column name or number of each sample's data frame that contains the measure of clonal abundance to be used when filtering the clusters. This column must have the same name in each sample.
}
  \item{min_seq_length}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{\link{buildRepSeqNetwork}} when performing
    network analysis on each sample (prior to filtering the clusters).
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{\link{buildRepSeqNetwork}}.
}
  \item{top_n_clusters}{
%%     ~~Describe \code{top_n_clusters} here~~
    The \code{top_n_clusters} clusters with the highest node count will
    automatically be included among the public clusters.
}
  \item{min_node_count}{
%%     ~~Describe \code{min_node_count} here~~
    Clusters with node count of at least this value will be included
    among the public clusters.
}
  \item{min_clone_count}{
%%     ~~Describe \code{min_clone_count} here~~
    Clusters with aggregate clone count of at least this value will be included
    among the public clusters.
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
    Passed to \code{\link{buildRepSeqNetwork}} when performing
    network analysis on each sample (prior to filtering the clusters).
}
  \item{print_plots}{
%%     ~~Describe \code{print_plots} here~~
    Passed to \code{\link{buildRepSeqNetwork}}.
}
  \item{plot_title}{
%%     ~~Describe \code{plot_title} here~~
    Passed to \code{\link{buildRepSeqNetwork}}.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Passed to \code{\link{buildRepSeqNetwork}}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    The outpout directory to save the filtered network results for each sample.
}
  \item{output_type}{
%%     ~~Describe \code{output_type} here~~
    A character string specifying the file format used to save the filtered
    results. Options include \code{"rds"}, \code{"csv"} and \code{"rda"}.
}
  \item{output_dir_unfiltered}{
%%     ~~Describe \code{output_dir_unfiltered} here~~
    The output directory to save the unfiltered network results for each sample,
    if desired. If \code{NULL} (the default), only the filtered results are
    saved.
}
  \item{output_type_unfiltered}{
%%     ~~Describe \code{output_type_unfiltered} here~~
    A character string specifying the file format scheme used to save the
    unfiltered network results for each sample. Only applicable if \code{output_dir_unfiltered}
    is non-null. Passed to \code{\link{buildRepSeqNetwork}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
    Other arguments to \code{buildRepSeqNetwork()}, not including \code{node_stats}, \code{stats_to_include} or \code{cluster_stats} (see details).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    All node-level and cluster-level network properties are automatically computed
    for each sample. The node-level properties
    are renamed to reflect their association to the sample-level network. Specifically,
    the properties are named
    \code{ClusterIDInSample}, \code{SampleLevelNetworkDegree},
    \code{SampleLevelTransitivity}, \code{SampleLevelCloseness},
    \code{SampleLevelCentralityByCloseness}, \code{SampleLevelCentralityByEigen},
    \code{SampleLevelEigenCentrality}, \code{SampleLevelBetweenness},
    \code{SampleLevelCentralityByBetweenness}, \code{SampleLevelAuthorityScore},
    \code{SampleLevelCoreness}, \code{SampleLevelPageRank}.
    A variable \code{SampleID} is
    added to both the node-level and cluster-level meta data for each sample.

    After the clusters in each sample are filtered, the node-level and cluster-level
    meta data are saved in the subdirectories \code{node_meta_data} and
    \code{cluster_meta_data}, respectively, of \code{output_dir}; each is saved
    using the sample ID (the corresponding element of \code{sample_id_list})
    as a filename prefix, followed by the appropriate file extension according
    to \code{output_type}.

    The files containing the node-level meta data for the filtered clusters can
    be supplied to \code{\link{buildPublicClusterNetwork}} in order to combine the
    public clusters across samples and perform network analysis. The files
    containing the cluster-level meta data for the filtered clusters can be
    supplied to \code{\link{buildPublicClusterNetworkByRepresentative}} to build a
    network using only a single representative sequence from each cluster.

    The unfiltered network results for each sample can also be saved by supplying
    a directory to \code{output_dir_unfiltered}, if these results are desired for
    downstream analysis.

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
\code{\link{buildPublicClusterNetwork}}

\code{\link{buildPublicClusterNetworkByRepresentative}}

\href{../doc/public_clusters.html}{Public Clusters Vignette}
}
\examples{
\dontrun{
vignette("public_clusters", package = "NAIR")
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
