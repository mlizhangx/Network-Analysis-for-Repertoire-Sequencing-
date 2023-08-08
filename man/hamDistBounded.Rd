\name{hamDistBounded}
\alias{hamDistBounded}
\title{
Bounded Computation of Hamming Distance
}
\description{
Computes the Hamming distance between two strings subject to a specified upper bound.
}
\usage{
hamDistBounded(a, b, k)
}
\arguments{
  \item{a}{
A character string.
}
  \item{b}{
A character string to be compared to \code{a}.
}
  \item{k}{
The upper bound on the Hamming distance between \code{a} and \code{b}.
}
}
\details{
For two character strings of equal length, the Hamming distance measures the total number of character differences between characters in corresponding positions. That is, for each position in one string, the character in that position is checked to see whether it differs from the character in the same position of the other string.

For two character strings of different lengths, the Hamming distance is not defined. However, the \code{hamDistBounded} function will accommodate strings of different lengths, doing so in a conservative fashion that seeks to yield a meaningful result for the purpose of checking whether two strings are sufficiently similar. If the two strings differ in length, placeholder characters are appended to the shorter string until its length matches that of the longer string. Each appended placeholder character is treated as different from the character in the corresponding position of the longer string. This is effectively the same as truncating the end of the longer string and adding the number of deleted characters to the Hamming distance between the shorter string and the truncated longer string (which is what is actually done in practice, as the computation is faster).

The above method used by \code{hamDistBounded} to accommodate unequal string lengths results in distance values whose meaning may be questionable, depending on context, when the two strings have different lengths. The decision to append placeholder characters to the end of the shorter string (as opposed to prepending them to the beginning) is ad hoc and somewhat arbitrary. In effect, it allows two strings of different lengths to be considered sufficiently similar if the content of the shorter string sufficiently matches the beginning content of the longer string and the difference in string length is not too great.

For comparing sequences of different lengths, the Levenshtein distance (see \code{\link{levDistBounded}}) is more appropriate and meaningful than using \code{hamDistBounded}, but comes at the cost of greater computational burden.

Computation is aborted early if the Hamming distance is determined to exceed the specified upper bound. This functionality is designed for cases when distinguishing between values above the upper bound is not meaningful, taking advantage of this fact to reduce the computational burden.
}
\value{
An integer. If the Hamming distance exceeds the specified upper bound \code{k}, then a value of \code{-1} is returned. Otherwise, returns the Hamming distance between \code{a} and \code{b}.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
\code{\link{levDistBounded}}
}
\examples{
# using an upper bound of 3
# (trivial since strings have length 3)
hamDistBounded("foo", "foo", 3)
hamDistBounded("foo", "fee", 3)
hamDistBounded("foo", "fie", 3)
hamDistBounded("foo", "foe", 3)
hamDistBounded("foo", "fum", 3)
hamDistBounded("foo", "bar", 3)

# using an upper bound of 1
# (most distances exceed the upper bound)
hamDistBounded("foo", "fee", 1)
hamDistBounded("foo", "fie", 1)
hamDistBounded("foo", "foe", 1)
hamDistBounded("foo", "fum", 1)
hamDistBounded("foo", "bar", 1)

# comparing strings of nonmatching length
hamDistBounded("foo", "fubar", 10)
hamDistBounded("foo", "foobar", 10)
hamDistBounded("foo", "barfoo", 10)

}