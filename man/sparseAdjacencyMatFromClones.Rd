\name{sparseAdjacencyMatFromClones}
\alias{sparseAdjacencyMatFromClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Sparse Adjacency Matrix for Immune Repertoire Network by Clone Sequence Similarity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a list of T-cell or B-cell receptor sequences, generates the adjacency matrix (in sparse format) for the corresponding network graph based on either the Hamming or Levenshtein distance using a specified cutoff.
}
\usage{
sparseAdjacencyMatFromClones(
    clones,
    dist_type = "hamming",
    max_dist = 1,
    drop_isolated_nodes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clones}{
%%     ~~Describe \code{clones} here~~
A character vector containing the TCR/BCR sequences.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
A character string specifying the distance type; valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{max_dist}{
%%     ~~Describe \code{max_dist} here~~
A positive integer specifying the maximum distance at which two TCR/BCR sequences are concidered adjacent (joined by an edge) in the network graph.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
A logical scalar; if \code{TRUE}, sequences with no adjacencies will have their rows/columns  excluded from the matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The adjacency matrix of a graph with \eqn{n} nodes/vertices is the symmetric \eqn{n \times n} matrix for which the \eqn{(i,j)}th entry is equal to 1 if nodes \eqn{i} and \eqn{j} are adjacent (joined by an edge in the graph) and 0 otherwise.

Each element of \code{clones} represents a node in the corresponding network graph; two nodes are adjacent if the similarity in their sequences, as measured by \code{dist_type}, is at most \code{max_dist}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A sparse matrix of class \code{dgCMatrix} from the \code{Matrix} package.

If \code{drop_isolated_nodes = TRUE}, the element of \code{clones} that corresponds to each row (equivalently, column) is specified in the \code{Dimnames} of this matrix. The \code{Dimnames} can be accessed using \code{dimnames()}. This returns a list with two elements, each a character vector whose length matches the row dimension of the matrix. These respectively specify the indices and the elements of \code{clones} corresponding to each row.
}
\references{
%% ~put references to the literature/web site here ~
    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{
sparseAdjacencyMatFromClones(
    c("fee", "fie", "foe", "fum", "foo"))
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
