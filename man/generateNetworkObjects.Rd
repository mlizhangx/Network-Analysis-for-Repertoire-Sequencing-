\name{generateNetworkObjects}
\alias{generateNetworkObjects}
\title{
    Generate Basic Output for an Immune Repertoire Network
}
\description{
    Given Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data, builds
    the network graph for the immune repertoire based on sequence similarity.
}
\usage{
generateNetworkObjects(data, seq_col,
                       dist_type = "hamming",
                       dist_cutoff = 1,
                       drop_isolated_nodes = TRUE)
}
\arguments{
  \item{data}{
    A data frame containing the AIRR-Seq data, with variables indexed by column
    and observations (e.g., clones or cells) indexed by row.
}
  \item{seq_col}{
    Specifies the column of \code{data} containing the receptor sequences to be
    used as the basis of similarity between rows. Accepts either the column name
    or column index. Alternatively, a vector of length 2 specifying two columns
    with distinct sequences (e.g., alpha chain and beta chain), in which case
    similarity between rows depends on both types of sequences being similar
    (e.g., similar alpha chain sequences and similar beta chain sequences).
}
  \item{dist_type}{
    Specifies the function used to measure the similarity between sequences.
    The similarity between two sequences determines the pairwise distance between
    their respective nodes in the network graph. Valid options are \code{"hamming"}
    (the default), which uses \code{\link{hamDistBounded}}, and \code{"levenshtein"},
    which uses \code{\link{levDistBounded}}. An additional option
    \code{"euclidean_on_atchley"} exists for TCR CDR3 amino acid sequences
    (see \code{\link{adjacencyMatAtchleyFromSeqs}}).
}
  \item{dist_cutoff}{
    A nonnegative scalar. Specifies the maximum pairwise distance (based on
    \code{dist_type}) for an edge connection to exist between two nodes. Pairs
    of nodes whose distance is less than or equal to this value will be joined
    by an edge connection in the network graph. Controls the stringency of the
    network construction and affects the number and density of edges in the network.
    A lower cutoff value requires greater similarity between sequences in order for
    their respective nodes to be joined by an edge connection. A value of \code{0}
    requires two sequences to be identical in order for their nodes to be joined
    by an edge.
}
  \item{drop_isolated_nodes}{
    A logical scalar. When \code{TRUE}, removes each node that is not joined by
    an edge connection to any other node in the network graph.
}
}
\details{
To construct the immune repertoire network, each TCR/BCR clone (bulk data)
or cell (single-cell data) is modeled as a node in the network graph,
corresponding to a single row of the AIRR-Seq data. For each node, the
corresponding receptor sequence is considered. Both nucleotide and amino acid
sequences are supported for this purpose. The receptor sequence is used as the
basis of similarity and distance between nodes in the network.

Similarity between sequences is measured using either the Hamming distance
or Levenshtein (edit) distance. The similarity determines the pairwise distance
between nodes in the network graph. The more similar two sequences are, the
shorter the distance between their respective nodes. Two nodes are joined by
an edge if their receptor sequences are sufficiently similar, i.e., if the
distance between the nodes is sufficiently small.

See the
\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html}{buildRepSeqNetwork}
package vignette for more details. The vignette can be accessed offline using
\code{vignette("buildRepSeqNetwork")}.
}
\value{
    If the constructed network contains no nodes, the function will return
    \code{NULL} with a warning. Otherwise, the function invisibly returns a list
    containing the following items:

    \item{igraph}{
    An object of class \code{\link[igraph]{igraph}} containing the list of nodes
    and edges for the network graph.
    }
    \item{adjacency_matrix}{
    The network graph adjacency matrix, stored as a sparse matrix of class
    \code{dgCMatrix} from the \code{Matrix} package. See
    \code{\link[Matrix]{dgCMatrix-class}}.
    }
    \item{node_data}{
    A data frame containing containing metadata for the network nodes, where each
    row corresponds to a node in the network graph. This data frame contains all
    variables from \code{data} (unless otherwise specified via \code{subset_cols})
    in addition to the computed node-level network properties if
    \code{node_stats = TRUE}. Each row's name is the name of the corresponding
    row from \code{data}.
    }
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html}{buildRepSeqNetwork vignette}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
toy_data <- simulateToyData()

net <-
  generateNetworkObjects(
    toy_data,
    "CloneSeq"
  )
}