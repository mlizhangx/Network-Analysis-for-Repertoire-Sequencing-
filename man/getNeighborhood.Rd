\name{getNeighborhood}
\alias{getNeighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Compute the Neighborhood Around a Cell/Clone
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given immune repertoire sequence data and a target receptor sequence, return
    the subset of the data consiting of those cells/clones whose receptor sequences are
    within a specified distance of the target.
}
\usage{
getNeighborhood(data, seq_col, target_seq,
                dist_type = "hamming", max_dist = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
    Data frame containing the RepSeq data.
}
  \item{seq_col}{
%%     ~~Describe \code{x} here~~
    The column name or number of \code{data} containing the receptor sequences of the same type as \code{target_seq}.
}
  \item{target_seq}{
%%     ~~Describe \code{x} here~~
    A character string containing the receptor sequence of the cell/clone around which to form a neighborhood.
}
  \item{dist_type}{
%%     ~~Describe \code{x} here~~
    A character string specifying the distance type on which the neighborhood is based; valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{max_dist}{
%%     ~~Describe \code{x} here~~
    A numeric scalar specifying the neighborhood radius. The neighborhood consists of all clones whose receptor sequences are at most \code{max_dist} from the target sequence.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The neighborhood will only be computed if \code{target_seq} is present in the data, i.e., in the column of \code{data} referenced by \code{seq_col}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If \code{target_seq} is present in the data, a data frame containing the subset corresponding to the neighborhood; otherwise, the function returns \code{NULL}.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Simulate some toy data for demonstration
toy_data <- simulateToyData(sample_size = 500)

# Get neighborhood around first clone sequence
nbd <- getNeighborhood(toy_data, seq_col = "CloneSeq",
                       target_seq = toy_data$CloneSeq[[1]])
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
