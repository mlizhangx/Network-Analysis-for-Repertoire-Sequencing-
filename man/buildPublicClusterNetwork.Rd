\name{buildPublicClusterNetwork}
\alias{buildPublicClusterNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Build Network of Public Clusters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given node-level meta data for each sample's filtered public clusters,
    combine the data across samples and perform network analysis.
}
\usage{
buildPublicClusterNetwork(

  ## Input ##
  file_list =
    list.files(file.path(getwd(), "public_clusters", "node_meta_data")),
  input_type = "rds", data_symbols = "ndat", header = TRUE, sep = "",
  seq_col,

  ## Network ##
  drop_isolated_nodes = FALSE,
  node_stats = NULL, stats_to_include = NULL, cluster_stats = NULL,

  ## Visualization ##
  color_nodes_by = "SampleID", color_scheme = "turbo",

  ## Output ##
  output_dir = file.path(getwd(), "public_clusters"),
  output_name = "PublicClusterNetwork",

  ...

)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} when loading the node-level meta data
  for each sample (one data frame per sample).
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
    The column name or number of each sample's data frame that contains the receptor sequences to be used as the basis of similarity between rows during network analysis. This column must have the same name in each sample.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{node_stats}{
%%     ~~Describe \code{node_stats} here~~
    Deprecated. Does nothing. All node-level network properties are automatically computed.
}
  \item{stats_to_include}{
%%     ~~Describe \code{stats_to_include} here~~
    Deprecated. Does nothing.
}
  \item{cluster_stats}{
%%     ~~Describe \code{cluster_stats} here~~
    Deprecated. Does nothing. All cluster-level network properties are automatically computed.
}
  \item{color_nodes_by}{
%%     ~~Describe \code{color_nodes_by} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{color_scheme}{
%%     ~~Describe \code{color_scheme} here~~
    Passed to \code{buildRepSeqNetwork}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
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

    The input data is intended to be obtained using \code{findPublicClusters()}
    on the unfiltered RepSeq data for each sample.

    Node-level network properties are renamed to reflect their association
    to the public network. The names used are
    \code{ClusterIDPublic}, \code{PublicNetworkDegree},
    \code{PublicTransitivity}, \code{PublicCloseness},
    \code{PublicCentralityByCloseness}, \code{PublicEigenCentrality},
    \code{PublicCentralityByEigen}, \code{PublicBetweenness},
    \code{PublicCentralityByBetweenness}, \code{PublicAuthorityScore},
    \code{PublicCoreness}, \code{PublicPageRank}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A list of network objects as returned by \code{buildRepSeqNetwork()}.
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
