\name{getNeighborhood}
\alias{getNeighborhood}
\title{
    Identify Cells or Clones in a Neighborhood Around a Target Sequence
}
\description{
    Given Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data and a
    target receptor sequence that is present within the data, identifies a
    "neighborhood" comprised of cells/clones with receptor sequences sufficiently
    similar to the target sequence.
}
\usage{
getNeighborhood(
    data,
    seq_col,
    target_seq,
    dist_type = "hamming",
    max_dist = 1
)
}
\arguments{
  \item{data}{
    A data frame containing the AIRR-Seq data.
}
  \item{seq_col}{
    Specifies the column of \code{data} containing the receptor sequences.
    Accepts either the column name or column index.
}
  \item{target_seq}{
    A character string containing the target receptor sequence. Must be a receptor
    sequence possessed by one of the clones/cells in the AIRR-Seq data.
}
  \item{dist_type}{
    Specifies the function used to quantify the similarity between receptor
    sequences. The similarity between two sequences determines their pairwise
    distance, with greater similarity corresponding to shorter distance. Valid
    options are \code{"hamming"} (the default), which uses
    \code{\link{hamDistBounded}}, and \code{"levenshtein"}, which uses
    \code{\link{levDistBounded}}.
}
  \item{max_dist}{
    Determines whether each cell/clone belongs to the neighborhood based on its
    receptor sequence's distance from the target sequence. The distance is based
    on the \code{dist_type} argument. \code{max_dist} specifies the maximum
    distance at which a cell/clone belongs to the neighborhood. Lower values
    require greater similarity between the target sequence and the receptor
    sequences of cells/clones in its neighborhood.
}
}
\value{
    A data frame containing the rows of \code{data} corresponding to the
    cells/clones in the neighborhood.

    If no cell/clone in the AIRR-Seq data possesses the target sequence as its
    receptor sequence, then a value of \code{NULL} is returned.
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

toy_data <- simulateToyData(sample_size = 500)

# Get neighborhood around first clone sequence
nbd <-
  getNeighborhood(
    toy_data,
    seq_col = "CloneSeq",
    target_seq = toy_data$CloneSeq[[1]]
  )

}