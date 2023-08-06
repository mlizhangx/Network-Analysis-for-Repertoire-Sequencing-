\name{aggregateIdenticalClones}
\alias{aggregateIdenticalClones}
\title{
    Aggregate Counts/Frequencies for Clones With Identical Receptor Sequences
}
\description{
    Given bulk Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) data with clones indexed by row, returns a data frame containing one row for each unique receptor sequence. Includes the number of clones sharing each sequence, as well as aggregate values for clone count and clone frequency across all clones sharing each sequence. Clones can be grouped according to metadata, in which case aggregation is performed within (but not across) groups.
}
\usage{
aggregateIdenticalClones(
  data,
  clone_col,
  count_col,
  freq_col,
  grouping_cols = NULL
)
}
\arguments{
  \item{data}{
      A data frame containing the bulk AIRR-Seq data, with clones indexed by row.
}
  \item{clone_col}{
      A column name or column index specifying the column of \code{data} containing the receptor sequences.
}
  \item{count_col}{
      A column name or column index specifying the column of \code{data} containing the clone counts.
}
  \item{freq_col}{
      A column name or column index specifying the column of \code{data} containing the clone frequencies.
}
  \item{grouping_cols}{
      An optional vector of column names or vector of column indices, specifying one or more columns of \code{data} used to assign clones to groups. If provided, aggregation occurs within groups, but not across groups. See details.
}
}
\details{
  If \code{grouping_cols} is left unspecified, the returned data frame will contain one row for each unique receptor sequence appearing in \code{data}.

  If one or more columns of \code{data} are specified using the \code{grouping_cols} argument, then each clone (row) in \code{data} is assigned to a group based on its combination of values in these columns. If two clones share the same receptor sequence but belong to different groups, their receptor sequence will appear multiple times in the returned data frame, with one row for each group in which the sequence appears. In each such row, the aggregate clone count, aggregate clone frequency, and number of clones sharing the sequence are reported within the group for that row.
}
\value{
    A data frame whose first column contains the receptor sequences and has the
    same name as the column of \code{data} specified by \code{clone_col}. One additional column will be present for each column of \code{data} that is specified using the \code{grouping_cols} argument, with each having the same column name. The remaining columns are as follows:

    \item{AggregatedCloneCount}{The aggregate clone count across all clones (within the same group, if applicable) that share the receptor sequence in that row.}

    \item{AggregatedCloneFrequency}{The aggregate clone frequency across all clones (within the same group, if applicable) that share the receptor sequence in that row.}

    \item{UniqueCloneCount}{The number of clones (rows) in \code{data} (within the same group, if applicable) possessing the receptor sequence for the current row.}
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \url{https://doi.org/10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
# Create some data
data <- data.frame(
  clone_seq = c("ATCG", rep("ACAC", 2), rep("GGGG", 4)),
  clone_count = rep(1, 7),
  clone_freq = rep(1/7, 7),
  time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
  subject_id = c(rep(1, 5), rep(2, 2))
)

# default usage
data_agg <-
  aggregateIdenticalClones(
    data, "clone_seq",
    "clone_count", "clone_freq"
  )

# group clones by time point
data_agg_time <-
  aggregateIdenticalClones(
    data, "clone_seq",
    "clone_count", "clone_freq",
    grouping_cols = "time_point"
  )

# group clones by subject ID
data_agg_subject <-
  aggregateIdenticalClones(
    data, "clone_seq",
    "clone_count", "clone_freq",
    grouping_cols = "subject_id"
  )

# group clones by time point and subject ID
data_agg_time_subject <-
  aggregateIdenticalClones(
    data, "clone_seq",
    "clone_count", "clone_freq",
    grouping_cols =
      c("subject_id", "time_point")
  )
}
