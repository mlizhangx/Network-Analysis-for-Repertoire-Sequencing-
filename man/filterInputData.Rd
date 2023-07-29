\name{filterInputData}
\alias{filterInputData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Filter RepSeq Data Rows and/or Subset Columns
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Filter data rows by sequence length and/or sequence content; subset data
    columns using a numeric or character vector of column references.
}
\usage{
filterInputData(data, seq_col, min_seq_length = NULL,
                drop_matches = NULL,
                subset_cols = NULL,
                count_col = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
    A data frame containing the immune repertoire sequencing data, with variables indexed by column and observations (e.g., clones or cells) indexed by row.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name of \code{data} containing the receptor sequences for use in any filtering based on length or content.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    A numeric scalar, or \code{NULL}. Observations whose receptor sequences have fewer than \code{min_seq_length} characters are removed prior to network analysis.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_chars} here~~
    An optional regular expression or character string. Observations whose receptor sequences return a match are removed prior to network analysis.
}
  \item{subset_cols}{
    Optional vector of column names or numbers of \code{data}; if supplied, only these (in addition to the column referenced by \code{seq_col}) are included in the output. If \code{NULL} (the default), all columns of \code{data} are included.
}
  \item{count_col}{
    Optional column name or number of \code{data}; if supplied, the specified column of the data will be coerced to numeric and data rows containing NA/NaNs in the count column will be dropped.
}
}
\details{
  The sequence column(s) of the data are coerced to character type and rows containing NA values in the sequence column(s) are dropped.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the filtered rows and specified columns of \code{data}.
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
## Generate some toy data for demonstration
raw_data <- simulateToyData()

# Filter data
filtered_data <-
  filterInputData(raw_data, seq_col = "CloneSeq",
                  min_seq_length = 13,
                  drop_matches = "[*|_]",
                  subset_cols = c("CloneSeq", "CloneFrequency", "SampleID"))
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