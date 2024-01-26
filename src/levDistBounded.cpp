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
int levDistBounded(std::string a,
                   std::string b,
                   const int& k) {

  if (k < 0) { return(-1); } // trivial bound
  if (a == b) { return(0); } // strings match
  if (k == 0) { return(-1); } // zero bound requires matching strings

  int n = a.length();
  int m = b.length();

  // if difference in string length exceeds bound, we are done
  if (abs(n - m) > k) { return(-1); }

  // If either string is empty, distance is the length of the other
  if (a.empty()) { return(m); }
  if (b.empty()) { return(n); }

  // Strip common prefix
  while (!a.empty() && !b.empty() && a.front() == b.front()) {
    a.erase(0, 1);
    b.erase(0, 1);
  }
  // Strip common suffix
  int bound = std::min(n, m), start;
  for (start = 0; start < bound && a[start] == b[start]; ++start);
  a = a.substr(start);
  b = b.substr(start);
  // Use shorter string for column dimension to save memory
  if (b.length() > a.length()) { a.swap(b); } // ensures b is not longer than a

  // Update string lengths
  n = a.length();
  m = b.length();

  if (b.empty()) { return(n); }

  // Bounded computation of Levenshtein distance:
  //   Full substring comparison matrix has |a|+1 rows and |b|+1 columns
  //   Compute only a diagonal band of width 2k+1, treating other values as Inf
  //   Only the current row and two scalars must be kept in memory
  int dist_inf = std::max(n, m) + 1; // effectively Inf
  Rcpp::IntegerVector current_row(m + 1);
  int prev_above;
  int prev_diag;
  int sub_cost;
  int ind_start;
  int ind_bound;
  bool bound_exceeded; // (distance bound; k)

  // fill row 0
  current_row.fill(dist_inf);
  int fill_bound = std::min(m, k) + 1;
  for (int i = 0; i < fill_bound; ++i) { current_row[i] = i; }
  // Rcout << "Row 0: \n" << current_row << "\n";

  // compute remaining rows
  for (int i = 1; i < (n + 1); ++i) {
    // Rcout << "======================================================\nProceeding to row " << i << ":\n";

    // indices of entries to be computed in diagonal band for current row
    ind_start = std::max(1, i - k); // column 0 already known
    ind_bound = std::min(m, i + k) + 1;
    // Rcout << "column indices of entries to be computed in current row: " << ind_start << " through " << stop_ind << "\n\n------------------------------------------\n\n";

    // Record value above-left of first band entry in new row,
    // i.e., first band entry of previous row.
    // This value will be assigned to prev_diag, with prev_above being updated,
    // before computing the first band entry in the new row.
    prev_above = current_row[ind_start - 1];
    // (values from previous row still stored in current_row)

    // If first band entry of current row is in the first column, it is
    // assigned the default value; if not, entry to its left treated as Inf
    if (i - k < 1) {
      current_row[0] = i;
    } else {
      current_row[ind_start - 1] = dist_inf;
    }

    // reset bound_exceeded
    bound_exceeded = true;

    // compute each band entry in current row
    for (int j = ind_start; j < ind_bound; ++j) {
      // Rcout << "computing entry in column " << j << ":\n";

      // does the current letter of a match current letter of b?
      sub_cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
      // Rcout << "Current letter of a: " << a[i - 1] << "\nCurrent letter of b: " << b[j - 1] << "\nAssociated substitution cost: " << sub_cost << "\n";

      // Update value above and to left of current entry
      prev_diag = prev_above;
      // Rcout << "Distance up to previous letter of a and previous letter of b (entry to upper-left): " << prev_diag << "\n";

      // Update value above current entry
      prev_above = current_row[j];
      // Rcout << "Distance up to previous letter of a (entry above): " << prev_above << "\n";

      // Rcout << "Distance up to previous letter of b (entry left): " << current_row[j - 1] << "\n";

      // update current entry
      current_row[j] = std::min(prev_above + 1,
                                std::min(current_row[j - 1] + 1,
                                         prev_diag + sub_cost));
      // Rcout << "Minimum distance selected; current entry updated to distance value of ** " << current_row[j] << " **\n\n";

      // if current entry is below bound, update bound_exceeded
      if (current_row[j] <= k) { bound_exceeded = false; }

    }

    // If all band entries in current row exceed bound, we are done
    if (bound_exceeded) { return(-1); }

  }

  // If final computed distance exceeds bound, return -1
  if (current_row[m] > k) { return(-1); }


  // return final computed distance if it is within bounds
  return(current_row[m]);

}