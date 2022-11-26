\name{getNeighborhood}
\alias{getNeighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Compute the Neighborhood Around a Cell/Clone
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given immune repertoire sequence data and a target receptor sequence, return
    the subset of the data consiting of those cells/clones whose receptor sequences are
    within a specified distance of the target.
}
\usage{
getNeighborhood(data, seq_col, target_seq,
                dist_type = "hamming", max_dist = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
    Data frame containing the RepSeq data.
}
  \item{seq_col}{
%%     ~~Describe \code{x} here~~
    The column name or number of \code{data} containing the receptor sequences of the same type as \code{target_seq}.
}
  \item{target_seq}{
%%     ~~Describe \code{x} here~~
    A character string containing the receptor sequence of the cell/clone around which to form a neighborhood.
}
  \item{dist_type}{
%%     ~~Describe \code{x} here~~
    A character string specifying the distance type on which the neighborhood is based; valid options are \code{"hamming"} and \code{"levenshtein"}.
}
  \item{max_dist}{
%%     ~~Describe \code{x} here~~
    A numeric scalar specifying the neighborhood radius. The neighborhood consists of all clones whose receptor sequences are at most \code{max_dist} from the target sequence.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The neighborhood will only be computed if \code{target_seq} is present in the data, i.e., in the column of \code{data} referenced by \code{seq_col}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If \code{target_seq} is present in the data, a data frame containing the subset corresponding to the neighborhood; otherwise, the function returns \code{NULL}.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Generate some data
set.seed(42)
sample_size <- 200
assign_group <- sample(c(TRUE, FALSE),
                       size = sample_size,
                       replace = TRUE)
base_seq_length <- 7
characters <- c("G", "A", "T", "C")
char_probs <- c(2/3, 2/9, 1/18, 1/18)
characters_sample <- sample(characters,
                            size = base_seq_length * sample_size,
                            replace = TRUE, prob = char_probs)
characters_sample <- matrix(characters_sample, ncol = sample_size)
clone_seqs <- apply(characters_sample, MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
clone_seqs_append <- rep("AATC", times = sample_size)
append_latent_prob <- runif(sample_size, min = 0, max = 1)
for (i in 1:sample_size) {
  case_group <- assign_group[[i]]
  if ((case_group & append_latent_prob[[i]] > 0.9) |
      (!case_group & append_latent_prob[[i]] > 0.5)) {
    clone_seqs_append[[i]] <- "AATCGGGG"
  } else if ((case_group & append_latent_prob[[i]] > 0.8) |
             (!case_group & append_latent_prob[[i]] > 0.3)) {
    clone_seqs_append[[i]] <- "AATCGGT"
  } else if ((case_group & append_latent_prob[[i]] > 0.7) |
             (!case_group & append_latent_prob[[i]] > 0.2)) {
    clone_seqs_append[[i]] <- "AATCGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.5) |
             (!case_group & append_latent_prob[[i]] > 0.1)) {
    clone_seqs_append[[i]] <- "AATTGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.25) |
             (!case_group & append_latent_prob[[i]] > 0.05)) {
    clone_seqs_append[[i]] <- "AATTG"
  }
  clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
                            collapse = "")
}
data <- data.frame(seqs = clone_seqs)
nbd <- getNeighborhood(data, seq_col = "seqs",
                       target_seq = data$seqs[[1]])
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
