\name{installPythonModules}
\alias{installPythonModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Install Python Modules Used By RepSeqNetworkAnalysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Checks whether the Python modules used by the `RepSeqNetworkAnalysis` package are available, and automatically installs any that are missing via `reticulate::py_install()`.
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
The modules checked for and installed are `numpy`, `scipy`, `pandas` and `levenshtein`. Calling this function is equivalent to making individual calls to `recitulate::py_install()` for each of these four modules, using the same arguments (except for `levenshtein`, which always uses `pip = TRUE`, as this module is not yet on Conda repositories).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Li Zhang, Hai Yang, Brian Neal
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
