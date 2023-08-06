\name{kmeansAtchley}
\alias{kmeansAtchley}
\title{
    Encode TCR Sequences Using Atchley Factors, Perform K-Means Clustering, Compare Samples by TCR Frequency Per Cluster
}
\description{
    Given multiple samples of AIRR-Seq data containing TCR CDR3 amino acid sequences, encodes each TCR sequence as a numeric vector based on Atchley factors and performs \eqn{K}-means clustering on the encoded sequences. Samples are compared based on the fraction of their unique TCR sequences in each cluster.
}
\usage{
kmeansAtchley(
  data,
  amino_col,
  sample_col,
  group_col,
  k = 100,
  pdf_width = 15,
  pdf_height = 15,
  margin_cluster_heatmap = 25,
  margin_corr_heatmap = 15,
  use_viridis = FALSE,
  output_dir = getwd(),
  file_cluster_heatmap =
    "atchley_kmeans_TCR_fraction_per_cluster.pdf",
  file_corr_heatmap =
    "atchley_kmeans_correlation_heatmap.pdf",
  return_output = FALSE
)
}
\arguments{
\item{data}{
    A data frame containing the combined AIRR-Seq data for all samples.
}
  \item{amino_col}{
    Specifies the column of \code{data} containing the TCR CDR3 amino acid sequence. Accepts either the column name or column index.
}
\item{sample_col}{
    Specifies the column of \code{data} containing the sample ID. Accepts either the column name or column index.
}
\item{group_col}{
    Specifies the column of \code{data} containing group IDs. Used to assign each sample to a group (e.g., treatment and control). Accepts either the column name or column index.
}
  \item{k}{
    The number of clusters to use for \eqn{K}-means clustering.
}
  \item{pdf_width}{
    Sets the page width when printing heatmaps to pdf. Passed to the \code{width} argument of the \code{\link[grDevices]{pdf}} function.
}
  \item{pdf_height}{
    Sets the page height when printing heatmaps to pdf. Passed to the \code{height} argument of the \code{\link[grDevices]{pdf}} function.
}
  \item{margin_cluster_heatmap}{
    The margin width of the column names and row names in the heatmap for each
    sample's fraction of TCRs per cluster.
}
  \item{margin_corr_heatmap}{
    The margin width of the column names and row names in the heatmap for the
    correlation in TCRs per cluster between samples.
}
  \item{use_viridis}{
    A logical scalar indicating whether to use the \code{\link[viridisLite]{cividis}} color scale instead of the \code{RdBu} color scale from \code{\link[RColorBrewer]{RColorBrewer}}. Better for greyscale printing and for viewers with common forms of color vision deficiency.
}
  \item{output_dir}{
    A file path. The directory in which to save pdf files of the two heatmaps produced.
}
  \item{file_cluster_heatmap}{
    A character string. The file name under which to save the pdf of the heatmap for each sample's fraction of TCRs per cluster.
}
  \item{file_corr_heatmap}{
    A character string. The file name under which to save the pdf of the heatmap for correlation in TCRs per cluster between samples.
}
  \item{return_output}{
    A logical scalar specifying whether to return the encoded TCR sequences, the K-means cluster membership table, and each sample's fraction of TCRs per cluster.
}
}
\details{
    Each unique TCR CDR3 amino acid sequence is encoded as a 30-dimensional numeric vector using \code{\link{encodeTCRSeqsByAtchleyFactor}}, which employs the software \href{https://doi.org/10.1038/s41592-020-01020-3}{TESSA}. \eqn{K}-means clustering is performed on the encoded TCR sequences from all samples, resulting in \code{k} clusters.

    For each sample, the frequency distribution of its TCR sequences across the \code{k} clusters is recorded in a length-\code{k} vector, where each entry is the fraction of TCRs in a particular cluster (with all \code{k} entries summing to 1). Fraction of TCRs refers to the fraction of the sample's total number of unique TCR sequences.

    The correlation in fraction of TCRs per cluster is computed between pairs of samples. This is done by pairing corresponding entries from both samples' \code{k}-vectors. Each of the \code{k} pairs of values is treated as an independent observation from a common bivariate random variable. This treats the fraction of TCRs in a cluster as a random variable, where the cluster is also considered random. It may be noted, however, that the \code{k} paired observations are not  actually independent, since they sum to \eqn{(1, 1)}.

    Two heatmaps are produced for comparing the samples. The first heatmap depicts each sample's fraction of TCRs in each cluster. The second heatmap depicts the pairwise correlation between samples that is described above. In both heatmaps, each column corresponds to a sample. A colored bar appears above each column, with the color indicating the group to which the sample belongs, based on the group assignments specified by the \code{group_ids} argument. This aids in the detection of differences and similarities within and between groups.
}
\value{
    If \code{return_output = TRUE}, then the function returns a list containing the following items:
    \item{kmeans_cluster_ids}{
    A data frame with two variables, \code{cdr3} and \code{kmeanClusterID}, containing the unique TCR sequences and the ID of the cluster to which each sequence belongs.
    }
    \item{encoded_values}{
    A matrix returned by \code{\link{encodeTCRSeqsByAtchleyFactor}} containing the numerically encoded vector for each TCR sequence.
    }
    \item{cluster_TCR_profiles}{
    A matrix with \code{k} rows (each corresponding to a cluster) and one column per sample. Each column contains a particular sample's fraction of TCRs in each cluster.
    }
}
\note{
The deep learning algorithm used by \href{https://doi.org/10.1038/s41592-020-01020-3}{TESSA} was trained specifically on TCR CDR3 amino acid sequences and is not appropriate for use with other amino acid sequences.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

    Zhang Z, Xiong D, Wang X, Liu H, Wang T. (2021). Mapping the functional landscape
of T cell receptor repertoires by single-T cell transcriptomics. \emph{Nat Methods}, 18
(1):92–9. \url{https://doi.org/10.1038/s41592-020-01020-3}

Atchley WR, Zhao J, Fernandes AD, Drüke T. (2005). Solving the protein sequence
metric problem. \emph{Proc Natl Acad Sci U.S.A.}  102(18):6395–400. \url{https://doi.org/10.1073/
pnas.0408677102}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
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