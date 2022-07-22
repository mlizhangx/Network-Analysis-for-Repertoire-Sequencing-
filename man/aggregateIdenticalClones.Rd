\name{aggregateIdenticalClones}
\alias{aggregateIdenticalClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Aggregate counts/frequencies for clones with identical receptor sequences
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
      A data frame containing the columns specified below.
}
  \item{clone_col}{
%%     ~~Describe \code{clone_col} here~~
      The column name or number of `data` containing the receptor sequences (e.g.,
      TCR CDR3 nucleotide sequence or amino acid sequence)
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
      The column name or number of `data` containing the clone counts
}
  \item{freq_col}{
%%     ~~Describe \code{freq_col} here~~
      The column name or number of `data` containing the clone frequencies
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
    A data frame whose first column has the same name as the column of `data` specified by `clone_col`
    containing the receptor sequences. If additional columns of `data` were supplied via the `grouping_cols`
    argument, these columns will also be inherited. The remaining columns are as follows:

    \item{AggregatedCloneCount}{The aggregated clone count for each receptor sequence}

    \item{AggregatedCloneFrequency}{The aggregated clone frequency for each receptor sequence}}
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

## The function is currently defined as
function (data, clone_col, count_col, freq_col, grouping_cols = NULL)
{
    if (is.numeric(clone_col)) {
        clone_col <- names(data)[clone_col]
    }
    if (is.numeric(count_col)) {
        count_col <- names(data)[count_col]
    }
    if (is.numeric(freq_col)) {
        freq_col <- names(data)[freq_col]
    }
    grouping_variables <- list(data[, clone_col])
    names(grouping_variables) <- clone_col
    if (!is.null(grouping_cols)) {
        if (is.numeric(grouping_cols)) {
            grouping_cols <- names(data)[grouping_cols]
        }
        for (i in 1:length(grouping_cols)) {
            grouping_variables$newvar <- data[, grouping_cols[[i]]]
            names(grouping_variables)[[length(grouping_variables)]] <- grouping_cols[[i]]
        }
    }
    cat("Aggregating reads (rows) by unique clone sequence...")
    data_to_aggregate <- list(AggregatedCloneCount = data[, c(count_col)],
        AggregatedCloneFrequency = data[, c(freq_col)])
    agg_counts <- stats::aggregate(data_to_aggregate, by = grouping_variables,
        FUN = sum)
    groups <- as.data.frame(grouping_variables)
    names(groups)[[1]] <- "temporary_placeholder_name"
    num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
        UniqueCloneCount = length(temporary_placeholder_name))
    names(num_reads)[[1]] <- clone_col
    out <- merge(agg_counts, num_reads, by = c(clone_col, grouping_cols))
    cat(paste0(" Done. ", nrow(out), " unique clone sequences found.\n"))
    return(out)
  }
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
