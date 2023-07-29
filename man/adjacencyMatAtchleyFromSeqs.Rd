\name{adjacencyMatAtchleyFromSeqs}
\alias{adjacencyMatAtchleyFromSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Network Adjacency Matrix Based on Atchley Factor Representations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
%%     ~~Describe \code{seqs} here~~
A character vector containing the TCR CDR3 amino acid sequences.
}
  \item{contig_ids}{
%%     ~~Describe \code{contig_ids} here~~
A numeric vector of the same length as \code{seqs}, used to uniquely identify its elements. By default, this is simply the sequence of integers from 1 to \code{length(seqs)}.
}
  \item{max_dist}{
%%     ~~Describe \code{max_dist} here~~
A positive integer specifying the maximum Euclidean distance at which the numeric vectors corresponding to two TCR sequences are concidered adjacent (joined by an edge) in the network graph.
}
  \item{return_type}{
%%     ~~Describe \code{return_type} here~~
A character string specifying the type of matrix to be returned. Valid options are \code{"adjacency_matrix"} and \code{"distance_matrix"}, the latter of which returns the matrix containing the pairwise Euclidean distances rather than the pairwise adjacencies.
}
  \item{outfile_distance_matrix}{
%%     ~~Describe \code{outfile_distance_matrix} here~~
An optional argument specifying a \code{csv} file to which the matrix will be written.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The adjacency matrix of a graph with \eqn{n} nodes/vertices is the symmetric \eqn{n \times n} matrix for which the \eqn{(i,j)}th entry is equal to 1 if nodes \eqn{i} and \eqn{j} are adjacent (joined by an edge in the graph) and 0 otherwise.

Each element of \code{seqs} represents a node in the corresponding network graph; two nodes are adjacent if the Euclidean distance between their numerically encoded values is at most \code{max_dist}.

If \code{return_type = "distance_matrix"}, then the \eqn{(i,j)}th entry of the matrix is the Euclidean distance between the numerically encoded vectors for the \eqn{i}th and \eqn{j}th sequences in the \code{seqs} vector.

The TCR CDR3 amino acid sequences are encoded as 30-dimensional numeric vectors using a trained encoder.
}
\value{
A symmetric, square matrix whose row/column dimension matches the length of \code{seqs}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\note{
The encoder was trained specifically on TCR CDR3 amino acid sequences and is not appropriate for use with other amino acid sequences.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
