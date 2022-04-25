# 0.0.9008
* New function `hamDistBounded` for computing bounded Hamming distance in C++
    * Supports strings of unequal length; the longer string is effectively truncated to the length of the shorter string, and the difference in length is added to the Hamming distance with the truncated version.  This is equivalent to extending the shorter string to the length of the longer string by appending placeholder characters that differ from their counterparts in the longer string.
* New function `hamAdjacencyMatSparse` for computing Hamming adjacency matrix in C++
* Changes to function `adjacencyMatrix`: 
    * Now supports `dist_type = "hamming"` to use Hamming distance for determining network adjacency
    * Original indices of input sequences that appear in the adjacency matrix are stored in the row names of the output matrix; the sequences themselves are stored in the column names (accessible via `dimnames()`)
    * The file `col_ids.txt` created by the C++ function that computes the adjacency matrix is now deleted by `adjacencyMatrix()` after it has finished its other tasks.  The information in the file is now stored in the row names of the output matrix and so the file is no longer needed.
    

# 0.0.9007
* `.genNetworkGraph()` internal helper for `buildNetwork()` renamed into a public version `genNetworkGraph()` for use by other package functions and by users; moved to a new file `utils.R` that will be used to house shared helper functions used by multiple package functions.
* Added `inst/python/Atchley_factors.csv`, which stores the Atchley factor amino acid embedding used by `BriseisEncoder.py`


# 0.0.9006
* Python module `tensorflow` added to `installPythonModules()` and Config/reticulate field of the DESCRIPTION file. `tensorflow` is required by the Python module `keras`.


# 0.0.9005
* file `zzz.R`  created with .onUnload() directive to unload package dll via call to `library.dynam.unload()` when the package is unloaded


# 0.0.9004
* Reintroduced reticulate package dependency as well as R .onLoad() directives and R functions related to python integration, dependency management and module installation
* Added Python source code for BriseisEncoder and h5 trained encoder file
* Removed $(SHLIB_OPENMP_CXXFLAGS) from PKG_CXXFLAGS and PKG_LIBS in Makevars as OpenMP is not supported for MacOS (still enabled in Makevars.win). Added comments with instructions for Linux users to enable OpenMP if desired, and added a corresponding note to the main Readme file's installation section for Linux users.


# 0.0.9003
* Removed Python source code (previously used to compute adjacency matrix)
* Removed reticulate package dependency (previously used to manage python dependencies)
* Removed R functions involved in Python management and calling Python scripts
* Changed internal helper R function .createAdjacencyMatrix() to public function adjacencyMatrix() and switched its calls to Python scripts for building the matrix into a call to the corresponding new C++ version added in 0.0.9002 (currently only Levenshtein distance is implemented; Hamming distance will be implemented in a future update)
* Changed buildNetwork() function to use adjacencyMatrix() instead of the previous function .createAdjacencyMatrix()
* Removed levAdjacencyMatDense in favor of always using levAdjacencyMatSparse
* Enabled ARMA_64BIT_WORD to accommodate larger matrices when using the C++ function levAdjacencyMatSparse; specifically, #define ARMA_64BIT_WORD=1 was added to the beginning of each .cpp source file, and -DARMA_64BIT_WORD=1 was added to the PKG_CXXFLAGS definition in Makevars and Makevars.win


# 0.0.9002
* Added internal support for Rcpp
* Added C++ routine for bounded Levenshtein distance (levDistBounded)
* Added C++ routine for dense Levenshtein adjacency matrix (levAdjacencyMatDense)
* Added C++ routine for sparse Levenshtein adjacency matrix (levAdjacencyMatSparse)


# 0.0.9001
Initial version
