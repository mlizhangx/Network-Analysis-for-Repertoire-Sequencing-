\name{installPythonModules}
\alias{installPythonModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    Install Python Modules Used by NAIR
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Checks whether the Python modules used by the `NAIR` package are available, and automatically installs any that are missing via `reticulate::py_install()`.
}
\usage{
  installPythonModules(method = "auto", conda = "auto", pip = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    %%     ~~Describe \code{method} here~~
      Passed to the `method` argument of `reticulate::py_install()`. See `?reticulate::py_install()` for details.
  }
  \item{conda}{
    %%     ~~Describe \code{conda} here~~
      Passed to the `conda` argument of `reticulate::py_install()`. See `?reticulate::py_install()` for details.
  }
  \item{pip}{
    %%     ~~Describe \code{pip} here~~
      Passed to the `pip` argument of `reticulate::py_install()`. See `?reticulate::py_install()` for details. For the `levenshtein` Python module, the call to `reticulate::py_install()` always uses `pip = TRUE`.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
    The modules checked for and installed are `numpy`, `pandas`, `tensorflow` and `keras`. Calling this function is equivalent to making individual calls to `recitulate::py_install()` for each of these modules.
}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
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
