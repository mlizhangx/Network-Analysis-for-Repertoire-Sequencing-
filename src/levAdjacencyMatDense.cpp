#include <RcppArmadillo.h>
#include <strings.h>
#include "levDistBounded.h"
using namespace arma;

// [[Rcpp::export]]
arma::umat levAdjacencyMatDense(std::vector<std::string> strings,
                                const int& maxdist) {

  // allocate memory for data structures
  const int dim = strings.size();
  int dist;
  arma::umat out = eye<umat>(dim, dim);  // identity matrix

  // compute adjacencies for upper triangle
  for (int j = 0; j < dim; ++j) {           // columns
    for (int i = 0; i < j ; ++i) {          // rows
      dist = levDistBounded(strings[i], strings[j], maxdist);
      if (dist != -1) { out(i, j) = 1; }
    }
  }

  // reflect upper triangle to lower
  out = arma::symmatu(out);

  // sum entries columnwise
  const arma::urowvec col_sums = arma::sum(out, 0);

  // record indices of nodes with positive degree
  const arma::uvec col_ids = find(col_sums > 1);

  // write indices of network nodes to file
  col_ids.save("col_ids.txt", raw_ascii);

  // subset matrix, keeping only nodes with positive degree
  out = out.submat(col_ids, col_ids);

  return(out);

}

