#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <strings.h>
#include "hamDistBounded.h"
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::sp_umat hamAdjacencyMatSparse(
    std::vector<std::string> strings,
    const int& maxdist,
    bool drop_deg_zero,
    std::string temp_dir
) {

  // allocate memory for data structures
  const int dim = strings.size();
  int dist;
  arma::sp_umat out = speye<sp_umat>(dim, dim);  // initialize as identity mat

  // compute adjacencies for upper triangle
  for (int j = 0; j < dim; ++j) {           // columns
    for (int i = 0; i < j; ++i) {           // rows
      dist = hamDistBounded(strings[i], strings[j], maxdist);
      if (dist != -1) { out(i, j) = 1; }
      Rcpp::checkUserInterrupt();
    }
  }

  // reflect upper triangle to lower
  out = arma::symmatu(out);

  if (drop_deg_zero) {

    // sum entries columnwise
    arma::sp_umat col_sums_spmat = arma::sum(out);
    arma::urowvec col_sums(col_sums_spmat);

    // record indices of nodes with positive degree
    arma::uvec col_ids = find(col_sums > 1);

    // subset matrix to keep only network nodes
    out = out.cols(col_ids);
    out = out.t();
    out = out.cols(col_ids);

    // write indices of network nodes to file
    col_ids += 1;  // offset C++'s 0-index starting convention
    col_ids.save(temp_dir + "/col_ids.txt", raw_ascii);

  }

  return(out);

}

