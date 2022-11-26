\name{adjacencyMatAtchleyFromSeqs}
\alias{adjacencyMatAtchleyFromSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Network Adjacency Matrix Based on Clone Sequence Atchley Factor Representations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a list of T-cell receptor CDR3 amino acid sequences, encodes each sequence as a 30-dimensional numeric vector based on the Atchley factor representation of its terms using a trained encoder, then computes the corresponding network adjacency matrix based on the Euclidean distance between encoded vectors using a specified cutoff.
}
\usage{
adjacencyMatAtchleyFromSeqs(
    seqs,
    contig_ids = seq_along(seqs),
    max_dist,
    return_type = "adjacency_matrix",
    outfile_distance_matrix = NULL)
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
}
\value{
A matrix.
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
%% ~put references to the literature/web site here ~
%%    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{
\dontrun{
cdr3 <- c("CASSEAQGSGSTDTQYF",
          "CATTEGSNTGELFF",
          "CASSIGDNEQFF",
          "CATSRDPDRGQSDTQYF",
          "CASSPTGLSGNTIYF",
          "CASSEEAGKDTQYF",
          "CASSGGADTQYF",
          "CASSLGLATDTQYF",
          "CASSEKEEVGELFF",
          "CASSSRTSGGAGELFF")
adjacencyMatAtchleyFromSeqs(
    cdr3, max_dist = 3)

adjacencyMatAtchleyFromSeqs(
    cdr3,
    return_type = "distance_matrix")
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
