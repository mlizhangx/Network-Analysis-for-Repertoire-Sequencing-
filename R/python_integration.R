# # Global references to python modules (values to be intialized by .onLoad)
# numpy <- pandas <- keras <- tensorflow <- NULL

# On Load:
.onLoad <- function(libname, pkgname) {

  # if a Python environment is already active, reticulate will configure the
  # active environment, including installing any packages in the
  # Config/reticulate field of the DESCRIPTION file
  reticulate::configure_environment(pkgname)

  # # Import python modules
  # numpy <<- reticulate::import("numpy", delay_load = TRUE)
  # pandas <<- reticulate::import("pandas", delay_load = TRUE)
  # keras <<- reticulate::import("keras", delay_load = TRUE)
  # tensorflow <<- reticulate::import("tensorflow", delay_load = TRUE)

  # ### DEPRECATED ###
  # # Offer to install miniconda and Python modules (deprecated; this is now handled automatically by reticulate)
  # request_msg <-
  #   "This package contains Python scripts and requires a distribution of Python 3.X to be installed in order to work properly. Would you like to install miniconda? (This requires 50-400 MB of storage space and may require time to download and install)."
  #
  # if (isTRUE(utils::askYesNo(request_msg, default = FALSE))) {
  #   reticulate::install_miniconda(force = TRUE)
  # }
  #
  #
  # message("NOTE: The following Python modules are required for this package to function properly: \n\n\t pandas \n\t scipy \n\t numpy \n\t Levenshtein \n\n You can install these using `RepSeqNetworkAnalysis::installPythonModules()`.\n")



}


# installs required Python modules
installPythonModules <- function(method = "auto", conda = "auto", pip = FALSE) {
  reticulate::py_install("numpy", method = method, conda = conda, pip = pip)
  reticulate::py_install("pandas", method = method, conda = conda, pip = pip)
  reticulate::py_install("tensorflow", method = method, conda = conda, pip = pip)
  reticulate::py_install("keras", method = method, conda = conda, pip = pip)
}



# Helper to check that required python modules are installed
.checkPythonModules <- function() {
  if (!reticulate::py_module_available("keras") |
      !reticulate::py_module_available("tensorflow") |
      !reticulate::py_module_available("numpy") |
      !reticulate::py_module_available("pandas")) {
    stop("One or more required Python modules are unavailable. Call `installPythonModules()` to install them.")
  }
}


# # Makes Python function for Briseis Encoder available in R
# initializeBriseisEncoder <- function() {
#   reticulate::source_python(system.file("python/BriseisEncoder_modified.py",
#                                         package = "RepSeqNetworkAnalysis")) }