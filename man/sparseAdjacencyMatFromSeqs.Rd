\name{sparseAdjacencyMatFromSeqs}
\alias{sparseAdjacencyMatFromSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute Network Adjacency Matrix for Immune Repertoire
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a list of receptor sequences, computes the adjacency matrix for the network graph based on sequence similarity.
}
\usage{
sparseAdjacencyMatFromSeqs(
    seqs,
    dist_type = "hamming",
    max_dist = 1,
    drop_isolated_nodes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
%%     ~~Describe \code{seqs} here~~
A character vector containing the receptor sequences.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
A character string specifying the distance function used to measure similarity between sequences. Valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{max_dist}{
%%     ~~Describe \code{max_dist} here~~
A nonnegative integer. Pairs of sequences whose distance is less than or equal to this value will be adjacent (connected by an edge) in the network graph.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
A logical scalar; if \code{TRUE}, sequences that are not adjacent to any other sequence will be excluded from the network.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The adjacency matrix of a graph with \eqn{n} nodes/vertices is the symmetric \eqn{n \times n} matrix for which the \eqn{(i,j)}th entry is equal to 1 if nodes \eqn{i} and \eqn{j} are adjacent (joined by an edge) in the network graph and 0 otherwise.

Each node in the network graph represents a receptor sequence from the \code{seqs} vector. Two nodes are adjacent if the distance between their sequences, as measured by \code{dist_type}, is less than or equal to \code{max_dist}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A sparse matrix of class \code{dgCMatrix} from the \code{Matrix} package.

If \code{drop_isolated_nodes = TRUE}, the row and column names of the matrix indicate which receptor sequences in the \code{seqs} vector correspond to each row and column. The row and column names can be accessed using \code{dimnames()}. This returns a list containing two character vectors, one for the row names and one for the column names. The name of the \eqn{i}th row is the numerical index of the \code{seqs} vector corresponding to the \eqn{i}th row and \eqn{i}th column. The  name of the \eqn{j}th column is the receptor sequence corresponding to the \eqn{j}th row and \eqn{j}th column.
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
\examples{
sparseAdjacencyMatFromSeqs(
    c("fee", "fie", "foe", "fum", "foo")
)
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
