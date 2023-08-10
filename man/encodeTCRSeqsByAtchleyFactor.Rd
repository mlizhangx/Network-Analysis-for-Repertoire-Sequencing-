\name{encodeTCRSeqsByAtchleyFactor}
\alias{encodeTCRSeqsByAtchleyFactor}
\title{
Encode TCR Sequences Numerically Based on Atchley Factor Representations
}
\description{
Given a list of T-cell receptor CDR3 amino acid sequences, encodes each sequence
as a 30-dimensional numeric vector based on the Atchley factor representations of
its amino acids.
}
\usage{
encodeTCRSeqsByAtchleyFactor(
    cdr3_AA,
    contig_ids = seq_along(cdr3_AA)
)
}
\arguments{
  \item{cdr3_AA}{
A character vector containing the TCR CDR3 amino acid sequences.
}
  \item{contig_ids}{
A numeric vector of the same length as \code{cdr3_AA}, used to uniquely identify
its elements. By default, this is simply the sequence of integers from 1 to
\code{length(cdr3_AA)}.
}
}
\details{
The encoding is performed using the software
\href{https://www.nature.com/articles/s41592-020-01020-3}{TESSA}. For each TCR
CDR3 amino acid sequence, each of its amino acids is represented by its
five-dimensional vector of Atchley factors, which numerically encode its biological
characteristics. The TCR CDR3 amino acid sequence is then encoded as a 30-dimensional
numeric vector based on the Atchley factor representations of its amino acids using
a deep learning algorithm.
}
\value{
A matrix with \code{length(cdr3_AA)} rows and 31 columns. Each row corresponds to
a TCR sequence; the first column is \code{contig_ids}, while the remaining columns
contain the 30-dimensional numerical encoding for each TCR sequence.
}
\note{
The deep learning algorithm used by
\href{https://www.nature.com/articles/s41592-020-01020-3}{TESSA} was trained
specifically on TCR CDR3 amino acid sequences and is not appropriate for use with
other amino acid sequences.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

Zhang Z, Xiong D, Wang X, Liu H, Wang T. (2021). Mapping the functional landscape
of T cell receptor repertoires by single-T cell transcriptomics.
\emph{Nat Methods}, 18(1):92–9.
\href{https://www.nature.com/articles/s41592-020-01020-3}{doi: 10.1038/s41592-020-01020-3}

Atchley WR, Zhao J, Fernandes AD, Druke T. (2005). Solving the protein sequence
metric problem. \emph{Proc Natl Acad Sci U.S.A.}  102(18):6395–400.
\href{https://www.pnas.org/doi/full/10.1073/pnas.0408677102}{doi: 10.1073/pnas.0408677102}


\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
\dontrun{
cdr3 <- c("CASSEAQGSGSTDTQYF",
          "CATTEGSNTGELFF",
          "CASSIGDNEQFF",
          "CATSRDPDRGQSDTQYF",
          "CASSPTGLSGNTIYF",
          "CASSEEAGKDTQYF",
          "CASSGGADTQYF",
          "CASSLGLATDTQYF",
          "CASSEKEEVGELFF",
          "CASSSRTSGGAGELFF")

encodeTCRSeqsByAtchleyFactor(cdr3)
}
}
