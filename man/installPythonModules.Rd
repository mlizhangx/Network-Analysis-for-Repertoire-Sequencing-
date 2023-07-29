\name{installPythonModules}
\alias{installPythonModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    Install Python Modules Used by NAIR
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Checks whether the Python modules used by the \code{NAIR} package are available, and automatically installs any that are missing via \code{reticulate::py_install()}.
}
\usage{
  installPythonModules(method = "auto", conda = "auto", pip = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    %%     ~~Describe \code{method} here~~
      Passed to the \code{method} argument of \code{reticulate::py_install()}. See \code{?reticulate::py_install()} for details.
  }
  \item{conda}{
    %%     ~~Describe \code{conda} here~~
      Passed to the \code{conda} argument of \code{reticulate::py_install()}. See \code{?reticulate::py_install()} for details.
  }
  \item{pip}{
    %%     ~~Describe \code{pip} here~~
      Passed to the \code{pip} argument of \code{reticulate::py_install()}. See \code{?reticulate::py_install()} for details. For the \code{levenshtein} Python module, the call to \code{reticulate::py_install()} always uses \code{pip = TRUE}.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
    The modules checked for and installed are \code{numpy}, \code{pandas}, \code{tensorflow} and \code{keras}. Calling this function is equivalent to making individual calls to \code{recitulate::py_install()} for each of these modules.
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


  % Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
  % \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
  % \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
