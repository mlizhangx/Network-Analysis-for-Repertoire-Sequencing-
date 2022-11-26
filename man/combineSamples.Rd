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
# Generate some data
set.seed(42)
num_samples <- 10; sample_size <- 10
groups <- sample(c(0, 1), size = num_samples, replace = TRUE)
init_chars <- c("AAAA", "AABA", "AACA", "AADA", "AAD", "AAA", "AABAA", "AAAAA")
init_probs_g0 <- rep(1 / length(init_chars), length(init_chars))
init_probs_g1 <- c(1, 5, 1, 1, 1, 1, 5, 1)
init_probs_g1 <- init_probs_g1 / sum(init_probs_g1)
samples <- matrix(nrow = sample_size, ncol = num_samples)
for (j in 1:num_samples) {
  init_probs <- init_probs_g0
  if (groups[[j]] == 1) { init_probs <- init_probs_g1 }
  samples[ , j] <- sample(init_chars, sample_size,
                          replace = TRUE, prob = init_probs)
}
num_edits <- 1
edit_ops <- c("insertion", "deletion", "transmutation")
edit_probs <- c(10/20, 2/20, 8/20)
chars <- c("A", "B", "C", "D")
char_probs_g0 <- c(1/2, 1/6, 1/6, 1/6)
char_probs_g1 <- c(1/3, 1/6, 1/6, 1/3)
for (k in 1:num_edits) {
  for (j in 1:num_samples) {
    char_probs <- char_probs_g0
    if ((groups[[j]] == 1)) { char_probs <- char_probs_g1 }
    for (i in 1:sample_size) {
      pos <- sample(1:nchar(samples[i, j]), size = 1)
      op <- sample(edit_ops, size = 1, prob = edit_probs)
      if (op == "insertion") {
        char <- sample(chars, size = 1, prob = char_probs)
        prefix <- substr(samples[i, j], 1, pos)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, char, suffix)
      } else if (op == "deletion") {
        prefix <- substr(samples[i, j], 1, pos - 1)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, suffix)
      } else if (op == "transmutation") {
        char <- sample(chars, size = 1, prob = char_probs)
        prefix <- substr(samples[i, j], 1, pos - 1)
        suffix <- substr(samples[i, j], pos + 1, nchar(samples[i, j]))
        samples[i, j] <- paste0(prefix, char, suffix)
      }
    }
  }
}
sample_symbols <- paste0("sample", 1:num_samples)
data_dir <- tempdir()
for (j in 1:num_samples) {
  dat <- data.frame(
    "seq" = samples[ , j], "sample_id" = j, "group" = groups[[j]],
    "count" = rbinom(sample_size, size = 100, prob = 0.5))
  assign(x = sample_symbols[[j]], value = dat)
  saveRDS(dat, file = file.path(data_dir, paste0(j, ".rds")))
  }


# Load data frames and combine
data <- loadDataFromFileList(
  file_list = file.path(data_dir, paste0(1:num_samples, ".rds")),
  input_type = "rds")


## Same as above, but filter out seqs shorter than 3 characters,
## drop the "count" column and add sample ID
data2 <- combineSamples(
  file_list = file.path(data_dir, paste0(1:num_samples, ".rds")),
  input_type = "rds",
  seq_col = "seq",
  min_seq_length = 3,
  subset_cols = "seq",
  sample_ids = as.character(1:num_samples))
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
