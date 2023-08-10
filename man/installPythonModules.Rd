\name{installPythonModules}
\alias{installPythonModules}
\title{
    Install Python Modules Used by NAIR
}
\description{
    Checks whether the Python modules used by the \code{NAIR} package are
    available, and automatically installs any that are missing via
    \code{\link[reticulate]{py_install}}.
}
\usage{
installPythonModules(
  method = "auto",
  conda = "auto",
  pip = FALSE
)
}
\arguments{
  \item{method}{
      Passed to the \code{method} argument of \code{\link[reticulate]{py_install}}.
  }
  \item{conda}{
      Passed to the \code{conda} argument of \code{\link[reticulate]{py_install}}.
  }
  \item{pip}{
      Passed to the \code{pip} argument of \code{\link[reticulate]{py_install}}.
  }
}
\details{
    The Python modules used by \code{NAIR} are \code{numpy}, \code{pandas},
    \code{tensorflow} and \code{keras}. Calling this function is equivalent to
    making individual calls to \code{\link[reticulate]{py_install}} for each of
    these modules.
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
