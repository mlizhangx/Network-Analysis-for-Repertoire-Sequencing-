\name{kmeansAtchley}
\alias{kmeansAtchley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Analyze TCR Samples Using Numerical Encoder and K-Means Clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given samples of TCR RepSeq data, convert the CDR3 amino acid
    sequences into numeric vectors using deep learning, perform \eqn{K}-means
    clustering, profile how each sample's unique TCR sequences are distributed
    among the clusters, and compare these profiles across samples.
}
\usage{
kmeansAtchley(
  data,
  amino_col = "AminoAcidSeq",
  sample_col = "SampleID",
  group_col, k = 100,
  pdf_width = 15,
  pdf_height = 15,
  margin_cluster_heatmap = 25,
  margin_corr_heatmap = 15,
  use_viridis = FALSE,
  output_dir = getwd(),
  file_cluster_heatmap =
    "atchley_kmeans_cluster_relative_size_profiles_by_sample.pdf",
  file_corr_heatmap =
    "atchley_kmeans_corr_in_cluster_size_profile_between_samples.pdf",
  return_output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
    A data frame containing the combined immune repertoire sequencing data for all samples,
    with variables indexed by column and observations indexed by row.
}
  \item{amino_col}{
%%     ~~Describe \code{amino_col} here~~
    The column of \code{data} containing the CDR3 amino acid sequences.
}
\item{sample_col}{
    The column of \code{data} containing the sample IDs.
}
\item{group_col}{
    The column of \code{data} containing the group IDs.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
    The number of clusters for \eqn{K}-means clustering.
}
  \item{pdf_width}{
%%     ~~Describe \code{plot_width} here~~
    Passed to the \code{width} argument of \code{grDevices::pdf()}.
}
  \item{pdf_height}{
%%     ~~Describe \code{plot_height} here~~
    Passed to the \code{height} argument of \code{grDevices::pdf()}.
}
  \item{margin_cluster_heatmap}{
%%     ~~Describe \code{margin_cluster_heatmap} here~~
    The margin width of column names and row names in the heatmap for each
    sample's cluster profile.
}
  \item{margin_corr_heatmap}{
%%     ~~Describe \code{margin_corr_heatmap} here~~
    The margin width of column names and row names in the heatmap for the
    correlation in cluster profiles between samples.
}
  \item{use_viridis}{
%%     ~~Describe \code{use_viridis} here~~
    A logical scalar indicating whether to use the \code{cividis}
    color scale from the \code{viridis} package instead of the \code{RdBu}
    color scale from \code{RColorBrewer}. This yields more readable results when
    printing in black-and-white, and is more robust against color
    vision deficiency.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    The directory in which to save pdf files of the two heatmaps produced.
}
  \item{file_cluster_heatmap}{
%%     ~~Describe \code{file_cluster_heatmap} here~~
    The pdf filename for the heatmap for each
    sample's cluster profile.
}
  \item{file_corr_heatmap}{
%%     ~~Describe \code{file_corr_heatmap} here~~
    The pdf filename for the heatmap for correlation in cluster profiles between samples.
}
  \item{return_output}{
%%     ~~Describe \code{return_output} here~~
    A logical scalar indicating whether to return the numeric vectors encoding
    the TCR sequences, the K-means cluster membership table, and each sample's
    profile of TCR representation across clusters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

    Each unique TCR sequence is encoded as a 30-dimensional numeric vector
    using \code{encodeTCRSeqsByAtchleyFactor()}, and \eqn{K}-means clustering is
    performed on the encoded TCR sequences from all samples, resulting in
    \code{k} clusters.

    For each sample, its TCR representation profile across the \code{k} clusters
    is computed as a length-\code{k} vector whose \eqn{i}th element is the
    fraction of the sample's unique TCR sequences that belong to the \eqn{i}th
    cluster.

    Two heatmaps are generated and saved to file: the first shows the TCR
    representation profile values across samples and clusters; the second shows
    the correlation in TCR represenation profiles between pairs of samples.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If \code{return_output = TRUE}, then a list containing the following items:
    \item{kmeans_cluster_ids}{
    A data frame with two variables, \code{cdr3} and \code{kmeanClusterID},
    containing the unique TCR sequences and the cluster to which each belongs.
    }
    \item{encoded_values}{
    A matrix returned by \code{encodeTCRSeqsByAtchleyFactor} containing the
    numerically encoded values of the TCR sequences.
    }
    \item{cluster_TCR_profiles}{
    A matrix with \code{k} rows and one column per sample; each column is the
    TCR representation profile for a sample, with each value recording the
    fraction of that sample's unique TCR sequences belonging to the cluster for
    the corresponding row (thus each column sums to 1).
    }
}
\references{
    The BriseisEncoder is authored by Ze Zhang (\email{Ze.Zhang@UTsouthwestern.edu}).
}

\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
\dontrun{
toy_data <- simulateToyData(
  samples = 20,
  sample_size = 50,
  prefix_length = 0,
  prefix_chars = "",
  prefix_probs = matrix(1, nrow = 20),
  affixes = c("CASSLGYEQYF", "CASSLGETQYF",
              "CASSLGTDTQYF", "CASSLGTEAFF",
              "CASSLGGTEAFF", "CAGLGGRDQETQYF",
              "CASSQETQYF", "CASSLTDTQYF",
              "CANYGYTF", "CANTGELFF",
              "CSANYGYTF"),
  affix_probs = matrix(1, ncol = 11, nrow = 20),
)
toy_data$GroupID <- rep(c("G1", "G2"), each = 500)

atchley <- kmeansAtchley(
  data = toy_data,
  k = 3,
  amino_col = "CloneSeq",
  sample_col = "SampleID",
  group_col = "GroupID",
  output_dir = tempdir(),
  return_output = TRUE)
}
}