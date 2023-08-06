\name{levDistBounded}
\alias{levDistBounded}
\title{
Bounded Computation of Levenshtein Distance
}
\description{
Computes the Levenshtein distance between two strings subject to a specified upper bound.
}
\usage{
levDistBounded(a, b, k)
}
\arguments{
  \item{a}{
A character string.
}
  \item{b}{
A character string to be compared to \code{a}.
}
  \item{k}{
The upper bound on the Levenshtein distance between \code{a} and \code{b}.
}
}
\details{
The Levenshtein distance (sometimes referred to as edit distance) between two character strings measures the minimum number of single-character edits (insertions, deletions and transformations) needed to transform one string into the other.

Compared to the Hamming distance (see \code{\link{hamDistBounded}}), the Levenshtein distance is particularly useful for comparing sequences of different lengths, as it can account for insertions and deletions, whereas the Hamming distance only accounts for single-character transformations. However, the computational burden for the Levenshtein distance can be significantly greater than for the Hamming distance.

Computation is aborted early if the Levenshtein distance is determined to exceed the specified upper bound. This functionality is designed for cases when distinguishing between values above the upper bound is not meaningful, taking advantage of this fact to reduce the computational burden.
}
\value{
An integer. If the Levenshtein distance exceeds the specified upper bound \code{k}, then a value of \code{-1} is returned. Otherwise, returns the Levenshtein distance between \code{a} and \code{b}.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{hamDistBounded}}
}
\examples{
# equal string lengths,
# character transmutations only
levDistBounded("foo", "bar", 3)
hamDistBounded("foo", "bar", 3) # agrees with Hamming distance

# one insertion, one deletion
levDistBounded("1234567", "1.23457", 7)
hamDistBounded("1234567", "1.23457", 7) # compare to Hamming distance

# same as above, but with a different lower bound
levDistBounded("1234567", "1.23457", 3) # within the bound
hamDistBounded("1234567", "1.23457", 3) # exceeds the bound

# one deletion (last position)
levDistBounded("1234567890", "123456789", 10)
hamDistBounded("1234567890", "123456789", 10)

# note the Hamming distance agrees with the Levenshtein distance
# for the above example, since the deletion occurs in the final
# character position. This is due to how hamDistBounded() handles
# strings of different lengths. In the example below, however...

# one deletion (first position)
levDistBounded("1234567890", "234567890", 10)
hamDistBounded("1234567890", "234567890", 10) # compare to Hamming distance

# one deletion, one transmutation
levDistBounded("foobar", "fubar", 6)
hamDistBounded("foobar", "fubar", 6) # compare to Hamming distance
}

