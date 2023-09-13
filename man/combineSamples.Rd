\name{combineSamples}
\alias{combineSamples}
\alias{loadDataFromFileList}
\title{
    Load and Combine Data From Multiple Samples
}
\description{
    Given multiple data frames stored in separate files,
    \code{loadDataFromFileList()}
    loads and combines them into a single data frame.

    \code{combineSamples()} has the same default behavior as
    \code{loadDataFromFileList()},
    but possesses additional arguments that allow the data frames to be filtered,
    subsetted and augmented with sample-level variables before being combined.
}
\usage{
loadDataFromFileList(
  file_list,
  input_type,
  data_symbols = NULL,
  header, sep, read.args
)

combineSamples(
  file_list,
  input_type,
  data_symbols = NULL,
  header, sep, read.args,
  seq_col = NULL,
  min_seq_length = NULL,
  drop_matches = NULL,
  subset_cols = NULL,
  sample_ids = NULL,
  subject_ids = NULL,
  group_ids = NULL,
  verbose = FALSE
)
}
\arguments{
    \item{file_list}{
    A character vector of file paths, or a list containing
    \code{\link[=connection]{connections}} and file paths.
    Each element corresponds to a single file containing the data
    for a single sample.
  }
    \item{input_type}{
    A character string specifying the file format of the sample data files.
    Options are \code{"rds"}, \code{"rda"}, \code{"csv"},
    \code{"csv2"}, \code{"tsv"}, \code{"table"}.
    See details.
  }
    \item{data_symbols}{
    Used when \code{input_type = "rda"}. Specifies the name of each sample's
    data frame within its respective Rdata file. Accepts a character vector of
    the same length as \code{file_list}. Alternatively, a single character string
    can be used if all data frames have the same name.
  }
  \item{header}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify a non-default value of the \code{header}
    argument to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
}
  \item{sep}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify a non-default value of the \code{sep}
    argument to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
}
  \item{read.args}{
    For values of \code{input_type} other than \code{"rds"} and \code{"rda"},
    this argument can be used to specify non-default values of optional
    arguments to \code{\link[utils:read.table]{read.table()}},
    \code{\link[utils:read.csv]{read.csv()}}, etc.
    Accepts a named list of argument values.
    Values of \code{header} and \code{sep}
    in this list take precedence over values specified via the \code{header}
    and \code{sep} arguments.
}
  \item{seq_col}{
    If provided, each sample's data will be filtered based on the values of
    \code{min_seq_length} and \code{drop_matches}.
    Passed to \code{\link[=filterInputData]{filterInputData()}} for each sample.
  }
    \item{min_seq_length}{
    Passed to \code{\link[=filterInputData]{filterInputData()}} for each sample.
  }
    \item{drop_matches}{
    Passed to \code{\link[=filterInputData]{filterInputData()}} for each sample.
  }
    \item{subset_cols}{
    Passed to \code{\link[=filterInputData]{filterInputData()}} for each sample.
  }
    \item{sample_ids}{
  A character or numeric vector of sample IDs, whose length matches that of
  \code{file_list}.
  }
    \item{subject_ids}{
  An optional character or numeric vector of subject IDs, whose length matches
  that of \code{file_list}. Used to assign a subject ID to each sample.
  }
    \item{group_ids}{
  A character or numeric vector of group IDs whose length matches that of
  \code{file_list}. Used to assign each sample to a group.
  }
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
}
\details{
    Each file is assumed to contain the data for a single sample, with
    observations indexed by row, and with the same columns across samples.

    Valid options for \code{input_type} (and the corresponding function used to
    load each file) include:

    \itemize{
    \item \code{"rds"}: \code{\link[=readRDS]{readRDS()}}
    \item \code{"rds"}: \code{\link[=readRDS]{readRDS()}}
    \item \code{"rda"}: \code{\link[=load]{load()}}
    \item \code{"csv"}: \code{\link[utils:read.csv]{read.csv()}}
    \item \code{"csv2"}: \code{\link[utils:read.csv]{read.csv2()}}
    \item \code{"tsv"}: \code{\link[utils:read.csv]{read.delim()}}
    \item \code{"table"}: \code{\link[utils:read.table]{read.table()}}
    }

    If \code{input_type = "rda"}, the \code{data_symbols} argument specifies the
    name of each data frame within its respective file.

    When calling \code{combineSamples()}, for each of \code{sample_ids},
    \code{subject_ids} and \code{group_ids} that is non-null, a corresponding
    variable will be added to the combined data frame; these variables are named
    \code{SampleID}, \code{SubjectID} and \code{GroupID}.
}
\value{
    A data frame containing the combined data rows from all files.
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
# Generate example data
set.seed(42)
samples <- simulateToyData(sample_size = 5)
sample_1 <- subset(samples, SampleID == "Sample1")
sample_2 <- subset(samples, SampleID == "Sample2")

# RDS format
rdsfiles <- tempfile(c("sample1", "sample2"), fileext = ".rds")
saveRDS(sample_1, rdsfiles[1])
saveRDS(sample_2, rdsfiles[2])

loadDataFromFileList(
  rdsfiles,
  input_type = "rds"
)

# With filtering and subsetting
combineSamples(
  rdsfiles,
  input_type = "rds",
  seq_col = "CloneSeq",
  min_seq_length = 13,
  drop_matches = "GGG",
  subset_cols = "CloneSeq",
  sample_ids = c("id01", "id02"),
  verbose = TRUE
)

# RData, different data frame names
rdafiles <- tempfile(c("sample1", "sample2"), fileext = ".rda")
save(sample_1, file = rdafiles[1])
save(sample_2, file = rdafiles[2])
loadDataFromFileList(
  rdafiles,
  input_type = "rda",
  data_symbols = c("sample_1", "sample_2")
)

# RData, same data frame names
df <- sample_1
save(df, file = rdafiles[1])
df <- sample_2
save(df, file = rdafiles[2])
loadDataFromFileList(
  rdafiles,
  input_type = "rda",
  data_symbols = "df"
)

# comma-separated values with header row; row names in first column
csvfiles <- tempfile(c("sample1", "sample2"), fileext = ".csv")
utils::write.csv(sample_1, csvfiles[1], row.names = TRUE)
utils::write.csv(sample_2, csvfiles[2], row.names = TRUE)
loadDataFromFileList(
  csvfiles,
  input_type = "csv",
  read.args = list(row.names = 1)
)

# semicolon-separated values with decimals as commas;
# header row, row names in first column
utils::write.csv2(sample_1, csvfiles[1], row.names = TRUE)
utils::write.csv2(sample_2, csvfiles[2], row.names = TRUE)
loadDataFromFileList(
  csvfiles,
  input_type = "csv2",
  read.args = list(row.names = 1)
)

# tab-separated values with header row and decimals as commas
tsvfiles <- tempfile(c("sample1", "sample2"), fileext = ".tsv")
utils::write.table(sample_1, tsvfiles[1], sep = "\t", dec = ",")
utils::write.table(sample_2, tsvfiles[2], sep = "\t", dec = ",")
loadDataFromFileList(
  tsvfiles,
  input_type = "tsv",
  header = TRUE,
  read.args = list(dec = ",")
)

# space-separated values with header row and NAs encoded as as "No Value"
txtfiles <- tempfile(c("sample1", "sample2"), fileext = ".txt")
utils::write.table(sample_1, txtfiles[1], na = "No Value")
utils::write.table(sample_2, txtfiles[2], na = "No Value")
loadDataFromFileList(
  txtfiles,
  input_type = "table",
  read.args = list(
    header = TRUE,
    na.strings = "No Value"
  )
)

# custom value separator and row names in first column
utils::write.table(sample_1, txtfiles[1],
                   sep = "@", row.names = TRUE, col.names = FALSE
)
utils::write.table(sample_2, txtfiles[2],
                   sep = "@", row.names = TRUE, col.names = FALSE
)
loadDataFromFileList(
  txtfiles,
  input_type = "table",
  sep = "@",
  read.args = list(
    row.names = 1,
    col.names = c("rownames",
                  "CloneSeq", "CloneFrequency",
                  "CloneCount", "SampleID"
    )
  )
)

# same as previous example
# (value of sep in read.args overrides value in sep argument)
loadDataFromFileList(
  txtfiles,
  input_type = "table",
  sep = "\t",
  read.args = list(
    sep = "@",
    row.names = 1,
    col.names = c("rownames",
                  "CloneSeq", "CloneFrequency",
                  "CloneCount", "SampleID"
    )
  )
)

\dontshow{
  # clean up temp directory
  unlink(c(rdsfiles, rdafiles, csvfiles, tsvfiles, txtfiles))
}
}