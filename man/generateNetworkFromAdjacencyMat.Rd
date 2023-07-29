\name{generateNetworkFromAdjacencyMat}
\alias{generateNetworkFromAdjacencyMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Generate the Network igraph for an Adjacency Matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  A simple wrapper function used for convenience. Given an adjacency matrix for
  an undirected network, generates the \code{igraph} object containing the undirected
  network graph corresponding to the adjacency matrix after simplifying and
  removing loops.
}
\usage{
generateNetworkFromAdjacencyMat(adjacency_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjacency_matrix}{
%%     ~~Describe \code{adjacency_matrix} here~~
    Passed to \code{igraph::graph_from_adjacency_matrix}. A square adjacency matrix.
    Sparse matrices created using the \code{Matrix} package are supported.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    An \code{igraph} graph object.
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
# Simulate some toy data for demonstration
toy_data <- simulateToyData(sample_size = 10)

adjmat <- sparseAdjacencyMatFromSeqs(toy_data$CloneSeq)
net <- generateNetworkFromAdjacencyMat(adjmat)

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
