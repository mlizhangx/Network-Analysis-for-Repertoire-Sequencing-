\name{generateNetworkObjects}
\alias{generateNetworkObjects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Generate Core Objects for a RepSeq Network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given immune repertoire sequencing data, returns the adjacency matrix,
    igraph and node meta data for the repetroire network based on sequence
    similarity.
}
\usage{
generateNetworkObjects(data, seq_col,
                       dist_type = "hamming",
                       dist_cutoff = 1,
                       drop_isolated_nodes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
    A data frame containing the RepSeq data, with clones/cells indexed by row.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of \code{data} containing the receptor sequences to be used as the basis of similarity between rows. Also accepts a vector of length 2 specifying distinct sequence columns (e.g., alpha, beta chains); then two rows are similar only if both types of sequences are similar.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    The type of function to use as a measure of similarity between two receptor
    sequences. Valid options are \code{"hamming"} (the default), \code{"levenshtein"}
    and \code{"euclidean_on_atchley"} (only applicable to TCR CDR3 amino acid sequences).
}
  \item{dist_cutoff}{
%%     ~~Describe \code{dist_cutoff} here~~
    A nonnegative scalar specifying the maximum distance threshold for
    similarity between receptor sequences.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
    A logical scalar; should observations whose receptor sequences are not
    similar to any other sequences be dropped from the network?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    This function is essentially a lightweight version of \code{buildRepSeqNetwork},
    omitting the input filtering, computation of network properties, generation
    of plots and saving of output.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If the constructed network contains no edges, the function will return \code{NULL} with a warning. Otherwise, a list containing the following items:

    \item{igraph}{An \code{igraph} object containing the edge list for the network.}
    \item{adjacency_matrix}{The network graph adjacency matrix, stored as a sparse matrix of class \code{dgCMatrix} from the \code{Matrix} package.}
    \item{node_data}{ A data
    frame containing the node-level meta-data for the network. This data frame contains
    all column names of \code{data}. The data frame will contain one row for
    each node that remains in the network after filtering based on
    \code{drop_isolated_nodes}, if applicable. The row names
    of the original input data will be preserved. }
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
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")
}
}