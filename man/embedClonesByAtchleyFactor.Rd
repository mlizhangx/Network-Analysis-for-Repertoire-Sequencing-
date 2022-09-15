\name{embedClonesByAtchleyFactor}
\alias{embedClonesByAtchleyFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Numerically Encode TCR CDR3 Amino Acid Sequences By Atchley Factor Representation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a list of T-cell receptor CDR3 amino acid sequences, encodes each sequence as a 30-dimensional numeric vector based on the Atchley factor representation of its terms using a trained encoder.
}
\usage{
embedClonesByAtchleyFactor(
    clones,
    contig_ids = seq_along(clones))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clones}{
%%     ~~Describe \code{clones} here~~
A character vector containing the TCR CDR3 amino acid sequences.
}
  \item{contig_ids}{
%%     ~~Describe \code{contig_ids} here~~
A numeric vector of the same length as \code{clones}, used to uniquely identify its elements. By default, this is simply the sequence of integers from 1 to \code{length(clones)}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Encoding is performed using the BriseisEncoder python script. Each amino acid is represented by its five-dimensional vector of Atchley factors, which numerically encode its biological characteristics. Each TCR CDR3 amino acid sequence is embedded in Euclidean 30-space based on the Atchley factor representations of its constitutent terms, using a deep learning encoder and a pretrained model.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix with \code{length(clones)} rows and 31 columns. Each row corresponds to a TCR sequence; the first column is \code{contig_ids}, while the remaining columns contain the 30-dimensional numerical encoding for each TCR sequence.
}
\note{
The encoder was trained specifically on TCR CDR3 amino acid sequences and is not appropriate for use with other amino acid sequences.
}
\references{
%% ~put references to the literature/web site here ~
    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{
\dontrun{
clones <- c("CASSEAQGSGSTDTQYF",
            "CATTEGSNTGELFF",
            "CASSIGDNEQFF",
            "CATSRDPDRGQSDTQYF",
            "CASSPTGLSGNTIYF",
            "CASSEEAGKDTQYF",
            "CASSGGADTQYF",
            "CASSLGLATDTQYF",
            "CASSEKEEVGELFF",
            "CASSSRTSGGAGELFF")
embedClonesByAtchleyFactor(clones)
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
