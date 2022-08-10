\name{hamDistBounded}
\alias{hamDistBounded}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Bounded Computation of Hamming Distance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes the Hamming distance between two strings subject to a specified upper bound.
}
\usage{
hamDistBounded(a, b, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
A character string.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
A character string to be compared to \code{a}.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
The upper bound on the Hamming distance between \code{a} and \code{b}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Computation is halted if the Hamming distance is determined to exceed the upper bound \code{k}; this reduces the computation necessary when distinguishing between values above the upper bound is unnecessary.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An integer whose value is the Hamming distance between \code{a} and \code{b} if this distance is at most \code{k}; otherwise the value is \code{-1}.
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
hamDistBounded("foo", "bar", 1)

hamDistBounded("foo", "bar", 10)
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
