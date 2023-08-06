


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