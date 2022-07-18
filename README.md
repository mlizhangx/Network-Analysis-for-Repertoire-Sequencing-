
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NAIR: Network Analysis for Immune Repeprtoire

<!-- badges: start -->
<!-- badges: end -->

This package implements pipelines for studying the immune repertoire of T and B cells 
via network analysis based on the similarity of their receptor sequences. 
NAIR can be used to relate clinical outcomes to immune repertoires based on their network 
properties, or to particular clusters and clones within a repertoire.

## Installation

You can install the development version of `NAIR` from
GitHub using the following commands in `R`:

install.packages("devtools")
devtools::install_github("mlizhangx/Network-Analysis-for-Repertoire-Sequencing-")

If you run into installation issues, consider updating `R` if using a
version prior to 3.0.2. Package installation requires
sufficiently-recent versions of the `R` packages `Rcpp` (1.0.8 or later)
and `RcppArmadillo` (10.8 or later); these version requirements are
automatically imposed via the `Imports` field of the `DESCRIPTION` file;
however, including minimum version requirements in the `Imports` field
is only supported in `R` from 3.0.2 onward.

**A note to Linux users:** If installing on Linux, note that OpenMP
support is disabled by default as MacOS does not support it. If you wish
to enable OpenMP, you can do so by editing the file src/Makevars;
editing instructions are included as comments in the file. OpenMP is
used by the Armadillo library for C++ to automatically parallelize
expensive operations such as elementwise matrix operations. Having
OpenMP enabled may shorten computation time for the network adjacency
matrix depending on the extent to which Armadillo’s methods for sparse
matrices take advantage of it.

### Python dependency

Certain functions within this package require Python version 3.6 or
greater. The dependency is handled automatically thanks to the
`reticulate` package; if a Python installation is not found when a
dependent function is called, the user will be automatically prompted to
install the latest miniconda Python distribution, and a Python
environment will be automatically set up.

`reticulate` will attempt to install any Python modules required by
functions in our package when setting up the Python environment. Package
functions that require Python modules to be available check first for
their availability. If they are unavailable, our package includes a
function `installPythonModules()` that can be called to check for and
dispatch calls to `reticulate::py_install()` in order to install any
missing modules.

If you wish to specify the Python environment to use, you can do so by
using the `use_virtualenv()` function from the `reticulate` package
prior to calling functions from `NAIR`.

## Vignettes

Pre-rendered vignettes are included in the package’s `vignettes`
directory; you can use the links below to view the rendered files on the
package’s GitHub repository. If you are installing the package from
GitHub wish to be able to access the vignettes from within `R`, include
the argument `build_vignettes = TRUE` when calling the
`devtools::install_github` function. Once the package is installed,
vignettes can be accessed from within `R` using
`help(package = "NAIR")` (to view in `R`) or
`browseVignettes(package = "NAIR")` (to view in web
browser).

[Generate Network With Node and Cluster
Statistics](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-/tree/main/vignettes/generateNetworkWithStats.md)
