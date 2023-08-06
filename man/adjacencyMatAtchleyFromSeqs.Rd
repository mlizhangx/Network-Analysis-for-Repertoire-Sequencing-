\name{adjacencyMatAtchleyFromSeqs}
\alias{adjacencyMatAtchleyFromSeqs}
\title{
Network Adjacency Matrix Based on Atchley Factor Representations
}
\description{
Given a list of T-cell receptor CDR3 amino acid sequences, encodes each sequence as a 30-dimensional numeric vector based on the Atchley factor representations of its amino acids, then computes the corresponding network adjacency matrix based on the Euclidean pairwise distances between encoded vectors.
}
\usage{
adjacencyMatAtchleyFromSeqs(
    seqs,
    contig_ids = seq_along(seqs),
    max_dist,
    return_type = "adjacency_matrix",
    outfile_distance_matrix = NULL
)
}
\arguments{
  \item{seqs}{
A character vector containing the TCR CDR3 amino acid sequences.
}
  \item{contig_ids}{
A numeric vector of the same length as \code{seqs}, used to uniquely identify its elements. By default, this is simply the sequence of integers from 1 to \code{length(seqs)}.
}
  \item{max_dist}{
A positive integer specifying the maximum Euclidean distance at which the numeric vectors corresponding to two TCR sequences are concidered adjacent (joined by an edge) in the network graph.
}
  \item{return_type}{
A character string specifying the type of matrix to be returned. Valid options are \code{"adjacency_matrix"} and \code{"distance_matrix"}, the latter of which returns the matrix containing the pairwise Euclidean distances rather than the pairwise adjacencies.
}
  \item{outfile_distance_matrix}{
An optional argument specifying a \code{csv} file to which the matrix will be written.
}
}
\details{
The adjacency matrix of a graph with \eqn{n} nodes/vertices is the symmetric \eqn{n \times n} matrix for which the \eqn{(i,j)}th entry is equal to 1 if nodes \eqn{i} and \eqn{j} are adjacent (joined by an edge in the graph) and 0 otherwise.

Each element of \code{seqs} represents a node in the corresponding network graph; two nodes are adjacent if the Euclidean distance between their numerically encoded values is at most \code{max_dist}.

If \code{return_type = "distance_matrix"}, then the \eqn{(i,j)}th entry of the matrix is the Euclidean distance between the numerically encoded vectors for the \eqn{i}th and \eqn{j}th sequences in the \code{seqs} vector.

The encoding is performed using \code{\link{encodeTCRSeqsByAtchleyFactor}}, which employs the software \href{https://doi.org/10.1038/s41592-020-01020-3}{TESSA}. For each TCR CDR3 amino acid sequence, each of its amino acids is represented by its five-dimensional vector of Atchley factors, which numerically encode its biological characteristics. The TCR CDR3 amino acid sequence is then encoded as a 30-dimensional numeric vector based on the Atchley factor representations of its amino acids using a deep learning algorithm.
}
\value{
A symmetric, square matrix (of class \code{Matrix}) whose row/column dimension matches the length of \code{seqs}.
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
# Example TCR CDR3 amino acid sequences
cdr3 <-
  c("CASSEAQGSGSTDTQYF",
    "CATTEGSNTGELFF",
    "CASSIGDNEQFF",
    "CATSRDPDRGQSDTQYF",
    "CASSPTGLSGNTIYF",
    "CASSEEAGKDTQYF",
    "CASSGGADTQYF",
    "CASSLGLATDTQYF",
    "CASSEKEEVGELFF",
    "CASSSRTSGGAGELFF"
  )

# Encode sequences,
# compute adjacency matrix
adjacencyMatAtchleyFromSeqs(
  cdr3,
  max_dist = 3
)

# Encode sequences,
# compute distance matrix
adjacencyMatAtchleyFromSeqs(
  cdr3,
  return_type = "distance_matrix"
)
}
}
