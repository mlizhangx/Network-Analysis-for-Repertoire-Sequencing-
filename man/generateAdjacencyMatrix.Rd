\name{generateAdjacencyMatrix}
\alias{generateAdjacencyMatrix}
\alias{sparseAdjacencyMatFromSeqs}
\title{
Compute Graph Adjacency Matrix for Immune Repertoire Network
}
\description{
Given a list of receptor sequences, computes the adjacency matrix for the
network graph based on sequence similarity.

\code{sparseAdjacencyMatFromSeqs()} is a deprecated equivalent of
\code{generateAdjacencyMatrix()}.
}
\usage{
generateAdjacencyMatrix(
  seqs,
  dist_type = "hamming",
  dist_cutoff = 1,
  drop_isolated_nodes = TRUE,
  verbose = FALSE
)

# Deprecated equivalent:
sparseAdjacencyMatFromSeqs(
  seqs,
  dist_type = "hamming",
  dist_cutoff = 1,
  drop_isolated_nodes = TRUE,
  verbose = FALSE,
  max_dist = deprecated()
)
}
\arguments{
  \item{seqs}{
A character vector containing the receptor sequences.
}
  \item{dist_type}{
Specifies the function used to quantify the similarity between sequences. The
similarity between two sequences determines the pairwise distance between their
respective nodes in the network graph, with greater similarity corresponding to
shorter distance. Valid options are \code{"hamming"} (the default), which uses
\code{\link{hamDistBounded}}, and \code{"levenshtein"}, which uses
\code{\link{levDistBounded}}.
}
  \item{dist_cutoff}{
A nonnegative scalar. Specifies the maximum pairwise distance (based on
\code{dist_type}) for an edge connection to exist between two nodes. Pairs of
nodes whose distance is less than or equal to this value will be joined by an
edge connection in the network graph. Controls the stringency of the network
construction and affects the number and density of edges in the network. A lower
cutoff value requires greater similarity between sequences in order for their
respective nodes to be joined by an edge connection. A value of \code{0} requires
two sequences to be identical in order for their nodes to be joined by an edge.
}
  \item{drop_isolated_nodes}{
Logical. When \code{TRUE}, removes each node that is not joined by an
edge connection to any other node in the network graph.
}
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link{stderr}}.
  }
  \item{max_dist}{
      \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
      Equivalent to \code{dist_cutoff}.
  }
}
\details{
The adjacency matrix of a graph with \eqn{n} nodes is the symmetric
\eqn{n \times n} matrix for which entry \eqn{(i,j)} is equal to 1 if nodes
\eqn{i} and \eqn{j} are connected by an edge in the network graph and 0 otherwise.

To construct the graph of the immune repertoire network, each receptor sequence
is modeled as a node. The similarity between receptor sequences, as measured
using either the Hamming or Levenshtein distance, determines the distance between
nodes in the network graph. The more similar two sequences are, the shorter the
distance between their respective nodes. Two nodes in the graph are joined by an
edge if the distance between them is sufficiently small, i.e., if their receptor
sequences are sufficiently similar.
}
\value{
A sparse matrix of class \code{dgCMatrix}
(see \code{\link[Matrix]{dgCMatrix-class}}).

If \code{drop_isolated_nodes = TRUE}, the row and column names of the matrix
indicate which receptor sequences in the \code{seqs} vector correspond to each
row and column of the matrix. The row and column names can be accessed using
\code{\link{dimnames}}. This returns a list containing two character vectors,
one for the row names and one for the column names. The name of the \eqn{i}th
matrix row is the index of the \code{seqs} vector corresponding to the \eqn{i}th
row and \eqn{i}th column of the matrix. The  name of the \eqn{j}th matrix column
is the receptor sequence corresponding to the \eqn{j}th row and \eqn{j}th column
of the matrix.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
generateAdjacencyMatrix(
  c("fee", "fie", "foe", "fum", "foo")
)

# No edge connections exist based on a Hamming distance of 1
# (returns a 0x0 sparse matrix)
generateAdjacencyMatrix(
  c("foo", "foobar", "fubar", "bar")
)

# Same as the above example, but keeping all nodes
# (returns a 4x4 sparse matrix)
generateAdjacencyMatrix(
  c("foo", "foobar", "fubar", "bar"),
  drop_isolated_nodes = FALSE
)

# Relaxing the edge criteria using a Hamming distance of 2
# (still results in no edge connections)
generateAdjacencyMatrix(
  c("foo", "foobar", "fubar", "bar"),
  dist_cutoff = 2
)

# Using a Levenshtein distance of 2, however,
# does result in edge connections
generateAdjacencyMatrix(
  c("foo", "foobar", "fubar", "bar"),
  dist_type = "levenshtein",
  dist_cutoff = 2
)

# Using a Hamming distance of 3
# also results in (different) edge connections
generateAdjacencyMatrix(
  c("foo", "foobar", "fubar", "bar"),
  dist_cutoff = 3
)
}