

# RepSeqNetworkAnalysis 0.0.9003

* Removed Python source code (previously used to compute adjacency matrix)
* Removed reticulate package dependency (previously used to manage python dependencies)
* Removed R functions involved in Python management and calling Python scripts
* Changed internal helper R function .createAdjacencyMatrix() to public function adjacencyMatrix() and switched its calls to Python scripts for building the matrix into calls to the corresponding new C++ versions added in 0.0.9002
* Changed buildNetwork() function to use adjacencyMatrix() instead of the previous function .createAdjacencyMatrix()
* Enabled ARMA_64BIT_WORD to accommodate larger matrices when using the C++ functions levAdjacencyMatDense and levAdjacencyMatSparse; specifically, #define ARMA_64BIT_WORD=1 was added to the beginning of each .cpp source file, and -DARMA_64BIT_WORD=1 was added to the PKG_CXXFLAGS definition in Makevars and Makevars.win


# RepSeqNetworkAnalysis 0.0.9002

* Added internal support for Rcpp
* Added C++ routine for bounded Levenshtein distance (levDistBounded)
* Added C++ routine for dense Levenshtein adjacency matrix (levAdjacencyMatDense)
* Added C++ routine for sparse Levenshtein adjacency matrix (levAdjacencyMatSparse)


# RepSeqNetworkAnalysis 0.0.9001

Initial version
