\name{filterInputData}
\alias{filterInputData}
\title{
    Filter Data Rows and Subset Data Columns
}
\description{
    Given a data frame with a column containing receptor sequences, filter data rows by sequence length and sequence content. Keep all data columns or choose which columns to keep.
}
\usage{
filterInputData(
    data,
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    count_col = NULL
)
}
\arguments{
  \item{data}{
    A data frame.
}
  \item{seq_col}{
    Specifies the column of \code{data} containing the receptor sequences. Accepts either the column name or column index. The column will be coerced to a character vector. Data rows containing a value of \code{NA} in this column will be dropped.
}
  \item{min_seq_length}{
    Observations whose receptor sequences have fewer than \code{min_seq_length} characters are dropped.
}
  \item{drop_matches}{
    Accepts a regular expression (see \code{\link[base]{regex}}). Checks values in the receptor sequence column for a pattern match using \code{\link{grep}}. Rows in which a match is found are dropped.
}
  \item{subset_cols}{
    Specifies which columns of the AIRR-Seq data are included in the output. Accepts a vector of column names or a vector of column indices. The default \code{NULL} includes all columns. The receptor sequence column is always included regardless of this argument's value.
}
  \item{count_col}{
    Specifies the column of \code{data} containing a variable with counts (such as clone count or UMI count). Accepts either the column name or column index. If provided, the specified column of the data will be coerced to numeric and data rows containing a value of \code{NA} or \code{NaN} in this column will be dropped.
}
}
\value{
    A data frame.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. \href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
raw_data <- simulateToyData()

# Remove sequences shorter than 13 characters,
# as well as sequences containing the subsequence "GGGG".
# Keep variables for clone sequence, clone frequency and sample ID
filtered_data <-
  filterInputData(
    raw_data,
    seq_col = "CloneSeq",
    min_seq_length = 13,
    drop_matches = "GGGG",
    subset_cols =
      c("CloneSeq", "CloneFrequency", "SampleID")
  )

}