
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    'Welcome to NAIR: Network Analysis of Immune Repertoire.\n',
    'Get started using `vignette("NAIR")`, or by visiting\n',
    'https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/'
  )
}