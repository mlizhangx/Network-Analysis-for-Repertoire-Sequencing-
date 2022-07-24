\name{filterClonesBySequenceLength}
\alias{filterClonesBySequenceLength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Filter Clones by Minimum Sequence Length
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given immune repertoire sequence data with clones indexed by row, remove clones
  with sequences below a specified minimum length
}
\usage{
filterClonesBySequenceLength(data, clone_col, min_length = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones
      indexed by row.
}
  \item{clone_col}{
%%     ~~Describe \code{clone_col} here~~
      The column name or number of `data` containing the sequence whose length
      is to be screened (e.g., TCR CDR3 nucleotide sequence or amino acid sequence).
      The column must be coercible to a character vector.
}
  \item{min_length}{
%%     ~~Describe \code{min_length} here~~
      A numeric value to be used as a cutoff; clones with sequences whose length
      is less than `min_length` will be removed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the same columns as `data` and those rows of `data`
    for which the sequence in `clone_col` has length at least `min_length`.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

  # Create some data
  data <- data.frame(
    clone_seq = c("ATCGATCG", rep("AC", 2), rep("GGGG", 4)),
    clone_count = rep(1, 7),
    clone_freq = rep(1 / 7, 7),
    # group/label variable 1
    time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
    # group/label variable 2
    subject_id = c(rep(1, 5), rep(2, 2))
  )

  # Filter using minimum sequence length 3 (default)
  data_filtered <-
    filterClonesBySequenceLength(data, "clone_seq")

  # Filter using minimum sequence length 5
  data_filtered5 <-
    filterClonesBySequenceLength(data, "clone_seq", min_length = 5)


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
