\name{combineSamples}
\alias{combineSamples}
\alias{loadDataFromFileList}
\title{
    Load and Combine Data From Multiple Samples
}
\description{
    Given multiple data frames stored in separate files, \code{loadDataFromFileList}
    loads and combines them into a single data frame.

    \code{combineSamples} is similar to \code{loadDataFromFileList},
    but allows the data frames to be filtered and subsetted before being
    combined, and can automatically add sample-level variables such as sample ID.
}
\usage{
loadDataFromFileList(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE, sep = ""
)

combineSamples(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE, sep = "",
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    sample_ids = NULL,
    subject_ids = NULL,
    group_ids = NULL
)
}
\arguments{
    \item{file_list}{
    A vector of file paths where each file contains the data for a sample, with
    observations indexed by row.
  }
    \item{input_type}{
    A character string specifying the file format of the sample data files.
    Options are \code{"table"}, \code{"txt"}, \code{"tsv"}, \code{"csv"},
    \code{"rds"} and \code{"rda"}.
  }
    \item{data_symbols}{
    Used when \code{input_type = "rda"}. Specifies the name of each sample's
    data frame within its respective Rdata file. Accepts a character vector of
    the same length as \code{file_list}. Alternatively, a single character string
    can be used if all data frames have the same name. See details.
  }
    \item{header}{
    Passed to \code{\link{read.table}} or \code{\link{read.csv}} for applicable
    values of \code{input_type}.
  }
    \item{sep}{
    Passed to \code{\link{read.table}} or \code{\link{read.csv}} for applicable
    values of \code{input_type}.
  }
  \item{seq_col}{
    Passed to \code{\link{filterInputData}} for each sample.
  }
    \item{min_seq_length}{
    Passed to \code{\link{filterInputData}} for each sample.
  }
    \item{drop_matches}{
    Passed to \code{\link{filterInputData}} for each sample.
  }
    \item{subset_cols}{
    Passed to \code{\link{filterInputData}} for each sample.
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
}
\details{
    Valid options for \code{input_type} (and the corresponding function used to
    load each file) include \code{"rds"} (\code{\link{readRDS}}), \code{"rda"}
    (\code{\link{load}}), \code{"csv"} (\code{\link{read.csv}}) and one of
    \code{"table"}, \code{"txt"} or \code{"tsv"} (all of which are loaded using
    \code{\link{read.table}}).

    If \code{input_type = "rda"}, the \code{data_symbols} argument specifies the
    name (in R) of each data frame as it is saved within its respective file.
    For example, given three R data frames named \code{sample1}, \code{sample2}
    and \code{sample3} in the R environment, with each data frame saved to a
    separate rda file using \code{\link{save}}), we would pass
    \code{c("sample1", "sample2", "sample3")} to \code{data_symbols}, assuming
    this is the order in which their file paths appear in the \code{file_list}
    argument. If each data frame is named \code{sample_data} within its
    respective rda file, we would simply specify
    \code{data_symbols = "sample_data"}.

    When calling \code{combineSamples}, for each of \code{sample_ids},
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
simulateToyData(sample_size = 5,
                output_dir = tempdir(),
                no_return = TRUE)

# Load and combine data
all_samples <-
  loadDataFromFileList(
    file_list =
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
    input_type = "rds",
  )

# Same as above, but also apply some filters,
# keep only the CloneSeq column and add SampleID column
all_samples <-
  combineSamples(
    file_list =
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
    input_type = "rds",
    seq_col = "CloneSeq",
    min_seq_length = 13,
    drop_matches = "GGG",
    subset_cols = "CloneSeq",
    sample_ids = c("id01", "id02")
  )


# clean up temp directory
file.remove(
  file.path(tempdir(), paste0("Sample", 1:2, ".rds"))
)
}