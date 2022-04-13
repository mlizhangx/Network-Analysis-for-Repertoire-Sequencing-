
# RepSeqNetworkAnalysis 0.0.9006
* Python module `tensorflow` added to `installPythonModules()` and Config/reticulate field of the DESCRIPTION file. `tensorflow` is required by the Python module `keras`.


# RepSeqNetworkAnalysis 0.0.9005
* file `zzz.R`  created with .onUnload() directive to unload package dll via call to `library.dynam.unload()` when the package is unloaded


# RepSeqNetworkAnalysis 0.0.9004

* Reintroduced reticulate package dependency as well as R .onLoad() directives and R functions related to python integration, dependency management and module installation
* Added Python source code for BriseisEncoder and h5 trained encoder file
* Removed $(SHLIB_OPENMP_CXXFLAGS) from PKG_CXXFLAGS and PKG_LIBS in Makevars as OpenMP is not supported for MacOS (still enabled in Makevars.win). Added comments with instructions for Linux users to enable OpenMP if desired, and added a corresponding note to the main Readme file's installation section for Linux users.


# RepSeqNetworkAnalysis 0.0.9003

* Removed Python source code (previously used to compute adjacency matrix)
* Removed reticulate package dependency (previously used to manage python dependencies)
* Removed R functions involved in Python management and calling Python scripts
* Changed internal helper R function .createAdjacencyMatrix() to public function adjacencyMatrix() and switched its calls to Python scripts for building the matrix into a call to the corresponding new C++ version added in 0.0.9002 (currently only Levenshtein distance is implemented; Hamming distance will be implemented in a future update)
* Changed buildNetwork() function to use adjacencyMatrix() instead of the previous function .createAdjacencyMatrix()
* Removed levAdjacencyMatDense in favor of always using levAdjacencyMatSparse
* Enabled ARMA_64BIT_WORD to accommodate larger matrices when using the C++ function levAdjacencyMatSparse; specifically, #define ARMA_64BIT_WORD=1 was added to the beginning of each .cpp source file, and -DARMA_64BIT_WORD=1 was added to the PKG_CXXFLAGS definition in Makevars and Makevars.win


# RepSeqNetworkAnalysis 0.0.9002

* Added internal support for Rcpp
* Added C++ routine for bounded Levenshtein distance (levDistBounded)
* Added C++ routine for dense Levenshtein adjacency matrix (levAdjacencyMatDense)
* Added C++ routine for sparse Levenshtein adjacency matrix (levAdjacencyMatSparse)


# RepSeqNetworkAnalysis 0.0.9001

Initial version
