\name{findAssociatedClones}
\alias{findAssociatedClones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Gather Clones in Public Neighborhoods Around Associated Sequences
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given multiple samples of RepSeq data from two groups and a list of receptor
    sequences associated to the comparison group (e.g., obtained using
    \code{findAssociatedSeqs}), this function does the following:
    for each associated sequence,
    from all samples that possess the associated sequence, gather the RepSeq data
    for all clones whose receptor sequences lie in a neighborhood of the associated sequence.

    The RepSeq data for each neighborhood's associated sequence is then stored in its own data frame
    and saved to its own file, in preparation for use with
    \code{buildAssociatedClusterNetwork}.
}
\usage{
findAssociatedClones(

  ## Input ##
  file_list, input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids = as.character(1:length(file_list)),
  subject_ids = NULL,
  group_ids,
  seq_col,

  ## Search Criteria ##
  assoc_seqs,
  nbd_radius = 1,
  dist_type = "hamming",
  min_seq_length = 6,
  drop_matches = "[*|_]",

  ## Output ##
  subset_cols = NULL,
  output_dir =
    file.path(getwd(), "associated_neighborhoods"),
  output_type = "csv",
  verbose = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_list}{
%%     ~~Describe \code{file_list} here~~
  Passed to \code{loadDataFromFileList} (one file per sample).
}
  \item{input_type}{
%%     ~~Describe \code{input_type} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{data_symbols}{
%%     ~~Describe \code{data_symbols} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{header}{
%%     ~~Describe \code{header} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
  Passed to \code{loadDataFromFileList}.
}
  \item{sample_ids}{
%%     ~~Describe \code{sample_ids} here~~
  A character or numeric vector of sample IDs, whose length matches that of \code{file_list}.
}
  \item{subject_ids}{
%%     ~~Describe \code{subject_ids} here~~
  An optional character or numeric vector of subject IDs, whose length matches that of \code{file_list}.
}
  \item{group_ids}{
%%     ~~Describe \code{group_ids} here~~
  A character or numeric vector of group IDs, which contains exactly two unique
  values and whose length matches that of \code{file_list}.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
    The column name or number of each sample's data frame that contains the receptor sequences to associate with the comparison group.
}
  \item{assoc_seqs}{
%%     ~~Describe \code{assoc_seqs} here~~
    A character vector containing the receptor sequences associated to the comparison group.
}
  \item{nbd_radius}{
%%     ~~Describe \code{nbd_radius} here~~
    The maximum distance between an element of \code{assoc_seqs} and other clones within its neighborhood.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    The distance metric used to measure the distance between sequences in \code{nbd_radius}. Valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{min_seq_length}{
%%     ~~Describe \code{min_seq_length} here~~
    Passed to \code{filterInputData()} when loading each sample.
}
  \item{drop_matches}{
%%     ~~Describe \code{drop_matches} here~~
    Passed to \code{filterInputData()}.
}
  \item{subset_cols}{
%%     ~~Describe \code{subset_cols} here~~
    Passed to \code{filterInputData()}.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
  The output directory; a valid output directory is required, since results are not returned in R.
}
  \item{output_type}{
%%     ~~Describe \code{output_dir} here~~
  A character string specifying the file format to use when writing each file. Valid options include \code{"csv"}, \code{"tsv"}, \code{"rds"} and \code{"rda"}.
}
  \item{verbose}{
%%     ~~Describe \code{output_dir} here~~
  A logical scalar; if \code{TRUE}, additional console output will be printed reporting the number of clones in each sample belonging to the neighborhood for each associated sequence, as well as the total number of clones in each neighborhood across all samples.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
   Each neighborhood's RepSeq data file is saved using the corresponding associated sequence
   for the filename, with the appropriate file extension appended based on the
   value of \code{output_type} (e.g., the neighborhood data for the sequence
   \code{"CASSGAYEQYF"} would be saved as \code{"CASSGAYEQYF.csv"} if
   \code{output_type = "csv"}).

   Three variables, \code{SampleID}, \code{SubjectID} and \code{GroupID}, are added to each neighborhood's RepSeq data, whose respective values in each row are the values of \code{sample_ids}, \code{subject_ids} and \code{group_ids} corresponding to the data in that row.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023). NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14. https://doi.org/10.3389/fimmu.2023.1181825

\url{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825}

\url{https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-}
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{findAssociatedSeqs}}
\code{\link{buildAssociatedClusterNetwork}}
}
\examples{
\dontrun{
## Generate some toy data for demonstration

# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_input_samples <- file.path(data_dir, "input_samples")
dir.create(dir_input_samples, showWarnings = FALSE)

samples <- 30
affixes <- c("AAAA", "AASA", "AACA", "AAQA", "AAQ", "AAA", "AASAA", "AAAAA")
affix_probs_g0 <- rep(1 / length(affixes),
                      times = length(affixes) * samples / 2)
affix_probs_g1 <- rep(c(1, 5, 1, 1, 1, 1, 5, 1), times = samples / 2)
affix_probs <- matrix(c(affix_probs_g0, affix_probs_g1),
                      nrow = samples, byrow = TRUE)
new_probs_g0 <- rep(c(1/2, 1/6, 1/6, 1/6), times = samples / 2)
new_probs_g1 <- rep(c(1/3, 1/6, 1/6, 1/3), times = samples / 2)
new_probs <- matrix(c(new_probs_g0, new_probs_g1),
                    nrow = samples, byrow = TRUE)
simulateToyData(
  samples = samples,
  sample_size = 30,
  prefix_length = 1,
  prefix_chars = c("A", "C"),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = affixes,
  affix_probs = affix_probs,
  num_edits = 4,
  edit_pos_probs = function(seq_length) {
    dnorm(seq(-4, 4, length.out = seq_length))
  },
  edit_ops = c("insertion", "deletion", "transmutation"),
  edit_probs = c(5, 1, 4),
  new_chars = c("A", "S", "C", "Q"),
  new_probs = new_probs,
  output_dir = dir_input_samples,
  no_return = TRUE
)


## 1. Find Associated Sequences

# input files for step 1 (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))
head(input_files)

# group labels for the samples
group_labels <- c(rep("reference", samples / 2), rep("comparison", samples / 2))

# search across samples for associated sequences using Fisher's exact test
associated_seqs <- findAssociatedSeqs(
  file_list = input_files, input_type = "rds",
  group_ids = group_labels, groups = c("reference", "comparison"),
  min_seq_length = NULL, drop_matches = NULL,
  seq_col = "CloneSeq", outfile = NULL)
head(associated_seqs)


## 2. Find Associated Clones

# output directory for current step
dir_nbds <- file.path(data_dir, "assoc_seq_nbds")

# Identify neighborhood around each associated sequence
findAssociatedClones(
  file_list = input_files, input_type = "rds", group_ids = group_labels,
  seq_col = "CloneSeq", dist_type = "levenshtein",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL, drop_matches = NULL,
  output_dir = dir_nbds)
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
