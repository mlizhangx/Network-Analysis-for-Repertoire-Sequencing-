\name{aggregateIdenticalClones}
\alias{aggregateIdenticalClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Aggregate Counts/Frequencies for Clones With Identical Receptor Sequences
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given TCR/BCR repertoire sequence data with clones indexed by row, aggregates
    rows with identical receptor sequences. Clone count and frequency are summed,
    other information is discarded and a column is added to record the number of
    unique clones (rows) for each unique receptor sequence.
}
\usage{
aggregateIdenticalClones(data, clone_col, count_col, freq_col,
                         grouping_cols = NULL)
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
      The column name or number of `data` containing the receptor sequences (e.g.,
      TCR CDR3 nucleotide sequence or amino acid sequence).
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
      The column name or number of `data` containing the clone counts.
}
  \item{freq_col}{
%%     ~~Describe \code{freq_col} here~~
      The column name or number of `data` containing the clone frequencies.
}
  \item{grouping_cols}{
%%     ~~Describe \code{grouping_cols} here~~
      An optional charcter or integer vector specifying one or more columns of `data`,
      each of which is treated as a group/label variable. If supplied, clones
      belonging to distinct groups will be treated as having distinct receptor
      sequences, effectively aggregating only the identical clones within each group.
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
    A data frame whose first column contains the receptor sequences and has the
    same name as the column of `data` specified by `clone_col`. If additional
    columns of `data` were supplied via the `grouping_cols` argument, these
    columns will also be inherited. The remaining columns are as follows:

    \item{AggregatedCloneCount}{The aggregated clone count for each receptor sequence.}

    \item{AggregatedCloneFrequency}{The aggregated clone frequency for each receptor sequence.}

    \item{UniqueCloneCount}{The number of clones (rows) in `data` possessing the receptor sequence for the current row
    (if groups are supplied, this is a within-group count).}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
  # Create some data
  data <- data.frame(
    clone_seq = c("ATCG", rep("ACAC", 2), rep("GGGG", 4)),
    clone_count = rep(1, 7),
    clone_freq = rep(1 / 7, 7),
    # group/label variable 1
    time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
    # group/label variable 2
    subject_id = c(rep(1, 5), rep(2, 2))
  )

  # Aggregate clones by receptor sequence (default usage)
  data_agg <- aggregateIdenticalClones(
    data, "clone_seq", "clone_count", "clone_freq")

  # Aggregate clones by receptor sequence and time point
  data_agg_time <- aggregateIdenticalClones(
    data, "clone_seq", "clone_count", "clone_freq",
    grouping_cols = "time_point")

  # Aggregate clones by receptor sequence and subject
  data_agg_subject <- aggregateIdenticalClones(
    data, "clone_seq", "clone_count", "clone_freq",
    grouping_cols = "subject_id")

  # Aggregate clones by receptor sequence, subject & time point
  # (note all clones in each group are already unique)
  data_agg_time_subject <- aggregateIdenticalClones(
    data, "clone_seq", "clone_count", "clone_freq",
    grouping_cols = c("subject_id", "time_point"))
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
