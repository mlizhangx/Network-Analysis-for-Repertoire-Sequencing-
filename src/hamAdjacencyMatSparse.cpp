// NAIR: Network Analysis of Immune Repertoire
// Copyright (C) 2023 Li Zhang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <strings.h>
#include "dropDegreeZero.h"
#include "hamDistBounded.h"
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".hamAdjacencyMatSparse")]]
arma::sp_umat hamAdjacencyMatSparse(
    std::vector<std::string> &strings,
    const int& maxdist,
    bool drop_deg_zero,
    std::string tempfile
) {
  // allocate memory for data structures
  std::unordered_map<std::string, std::vector<int>> str2idx; // keep map from unique strings to indices
  const int dim = strings.size();
  int dist;
  arma::sp_umat out = speye<sp_umat>(dim, dim);  // initialize as identity mat

  for (int i = 0; i < dim; i++)
    str2idx[strings[i]].push_back(i);

  // compute adjacencies
  for (auto it1 = str2idx.begin(); it1 != str2idx.end(); ++it1) {
    for (auto it2 = str2idx.begin(); it2 != it1; ++it2) {
      dist = hamDistBounded(it1->first, it2->first, maxdist);
      if (dist != -1)
        for (auto idx1 : it1->second)
          for (auto idx2 : it2->second)
            out(idx1, idx2) = out(idx2, idx1) = 1; 
    }
    for (auto idx1 : it1->second)
      for (auto idx2 : it1->second)
        out(idx1, idx2) = out(idx2, idx1) = 1; 
    Rcpp::checkUserInterrupt();
  }

  dropDegreeZero(drop_deg_zero, out, tempfile);
  return(out);

}

