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
#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;

// [[Rcpp::export]]
int hamDistBounded(const std::string &a,
                   const std::string &b,
                   const int& k) {

  if (k < 0) { return(-1); } // trivial bound
  if (a == b) { return(0); } // strings match
  if (k == 0) { return(-1); } // zero bound requires matching strings

  int n = a.length();
  int m = b.length();

  // Initialize distance value, bound below by difference in string length
  int dist = abs(n - m);

  // stop if distance exceeds bound
  if (dist > k) { return(-1); }

  // Compute hamming distance; longer string truncated to length of shorter
  int ind_bound = std::min(n, m);
  for (int i = 0; i < ind_bound; ++i) {
    if (a[i] != b[i]) { 
      dist++; 
      if (dist > k) 
        return(-1); // stop if distance exceeds bound
    }
  }

  // return distance
  return(dist);

}