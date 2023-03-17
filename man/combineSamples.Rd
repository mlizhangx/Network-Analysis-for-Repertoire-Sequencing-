\name{combineSamples}
\alias{combineSamples}
\alias{loadDataFromFileList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Load Multiple Data Files and Combine
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple data frames stored in separate files, \code{loadDataFromFileList}
    loads and combines them into a single data frame.

    \code{combineSamples} is similar to \code{loadDataFromFileList},
    but allows the data frames to be filtered and subsetted before being
    combined, and can automatically add sample-level variables such as sample ID.
}
\usage{
loadDataFromFileList(
    file_list, input_type,
    data_symbols = NULL,
    header = TRUE, sep = "")

combineSamples(
    file_list, input_type,
    data_symbols = NULL,
    header = TRUE, sep = "",
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    sample_ids = NULL,
    subject_ids = NULL,
    group_ids = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{file_list}{
  %%     ~~Describe \code{file_list} here~~
    A character vector of file paths containing the data frames (one file per data frame).
  }
    \item{input_type}{
  %%     ~~Describe \code{input_type} here~~
    A character string specifying the input file format; determines the function used to load each file. See details section below.
  }
    \item{data_symbols}{
  %%     ~~Describe \code{data_symbols} here~~
    If \code{input_type = "rda"}, a character vector specifying the name of each data frame (this can be of length 1 if all data frames have the same name).
  }
    \item{header}{
  %%     ~~Describe \code{header} here~~
    Passed to \code{read.table} or \code{read.csv} if applicable.
  }
    \item{sep}{
  %%     ~~Describe \code{sep} here~~
    Passed to \code{read.table} or \code{read.csv} if applicable.
  }
  \item{seq_col}{
  %%     ~~Describe \code{seq_col} here~~
    Passed to \code{filterInputData} for each sample.
  }
    \item{min_seq_length}{
  %%     ~~Describe \code{min_seq_length} here~~
    Passed to \code{filterInputData} for each sample.
  }
    \item{drop_matches}{
  %%     ~~Describe \code{drop_matches} here~~
    Passed to \code{filterInputData} for each sample.
  }
    \item{subset_cols}{
  %%     ~~Describe \code{subset_cols} here~~
    Passed to \code{filterInputData} for each sample.
  }
    \item{sample_ids}{
  %%     ~~Describe \code{sample_ids} here~~
  An optional character or numeric vector of sample IDs, whose length matches that of \code{file_list}.
  }
    \item{subject_ids}{
  %%     ~~Describe \code{subject_ids} here~~
  An optional character or numeric vector of subject IDs, whose length matches that of \code{file_list}.
  }
    \item{group_ids}{
  %%     ~~Describe \code{group_ids} here~~
  An optional character or numeric vector of group IDs, whose length matches that of \code{file_list}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    Valid options for \code{input_type} (and the corresponding function used to
    load each file) include \code{"rds"} (\code{readRDS}), \code{"rda"}
    (\code{load}), \code{"csv"} (\code{read.csv}) and \code{"table"}
    (\code{read.table}).

    If \code{input_type = "rda"}, the \code{data_symbols}
    argument specifies the name of each data frame in R. For example,
    given three R data frames named \code{sample1}, \code{sample2}
    and \code{sample3}, with each data frame saved to a separate rda file using
    \code{save()}), we would pass \code{c("sample1", "sample2", "sample3")}
    to \code{data_symbols}.

    For \code{combineSamples}, for each of \code{sample_ids}, \code{subject_ids}
    and \code{group_ids} that is non-null, a corresponding variable will be
    added to the combined data frame; these will be named \code{SampleID},
    \code{SubjectID} and \code{GroupID}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing the combined data rows from all files.
}
\references{
%% ~put references to the literature/web site here ~
%%    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{
## Generate some toy data for demonstration

# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_samples <- file.path(data_dir, "samples")
dir.create(dir_samples, showWarnings = FALSE)
simulateToyData(sample_size = 5,
                output_dir = dir_samples,
                no_return = TRUE)


# Load data frames and combine
data <- loadDataFromFileList(
  file_list = list.files(dir_samples, full.names = TRUE),
  input_type = "rds")


## Same as above, but filter out seqs shorter than 3 characters,
## drop the count/freq/sampleID columns & add new sample ID column
data2 <- combineSamples(
  file_list = list.files(dir_samples, full.names = TRUE),
  input_type = "rds",
  seq_col = "CloneSeq",
  min_seq_length = 3,
  subset_cols = "CloneSeq",
  sample_ids = c("id01", "id02"))
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
