//compute Levenshtein distance between two strings a and b of lengths n and m,
//subject to a specified upper bound k.  The upper bound allows faster
//computation since, rather than dynamically computing the entire (n+1 x m+1)
//matrix of Levenshtein distances on partial strings of a and b, we can restrict
//attention to a diagonal band whose width depends on k.  If it is found that
//the Levenshtein distance will exceed the bound, the computation is halted and
//a value of -1 is returned.

#include <Rcpp.h>
#include <strings.h>
using namespace Rcpp;

// [[Rcpp::export]]
int levDistBounded(std::string a,
                   std::string b,
                   const int& k) {

  // if strings are identical, return 0
  if (a == b) { return(0); }

  // get string lengths
  int n = a.length();
  int m = b.length();
  // Rcout << "length of a: " << n << " characters\n";
  // Rcout << "length of b: " << m << " characters\n";

  // lev dist is bounded below by difference in string length;
  // if this difference exceeds the specified bound, we are done
  if (abs(n - m) > k) { return(-1); }

  // Strip common prefix
  while (a[0] == b[0]) {
    // Rcout << "First character of both strings match: Removing character " << a[0] << "...\n";
    a = a.substr(1, n - 1);
    b = b.substr(1, m - 1);
    n = n - 1;
    m = m - 1;
    // Rcout << "a is now " << a << "\n";
    // Rcout << "b is now " << b << "\n";
  }

  // Strip common suffix
  while (a[n - 1] == b[m - 1]) {
    // Rcout << "Last character of both strings match: Removing character " << a[a.length() - 1] << "...\n";
    a = a.substr(0, n - 1);
    b = b.substr(0, m - 1);
    n = n - 1;
    m = m - 1;
    // Rcout << "a is now " << a << "\n";
    // Rcout << "b is now " << b << "\n";
  }

  // Use shorter string for row dimension to save memory
  if (m > n) {
    // Rcout << "b is longer than a; swapping to save memory...\n";
    std::string tmp = a;
    a = b;
    b = tmp;
    int tmp2 = m;
    m = n;
    n = tmp2;
  }

  // if b is empty, the distance is the length of a
  if (m == 0) { return(n); }

  // Bounded computation of Levenshtein distance:
  //   compute diagonal band of matrix entries
  //   band extends k entries to either side of each main diagonal entry
  //   treat entries outside of this band as infinite-valued

  // Get value to use as effective "infinite distance"
  int dist_inf = std::max(n, m) + 1;

  // allocate memory for current row
  Rcpp::IntegerVector current_row(m + 1);

  // allocate memory for two entries of previous row and other values
  int prev_above;
  int prev_diag;
  int sub_cost;
  int start_ind;
  int stop_ind;
  bool bound_exceeded;

  // populate current row with values of first row
  current_row.fill(dist_inf);
  const int fill_bound = 1 + std::min(m, k);
  for (int i = 0; i < fill_bound; ++i) current_row[i] = i;

  // Rcout << "Row 0: \n" << current_row << "\n";

  // compute values in diagonal band for each row beyond the first
  for (int i = 1; i < (n + 1); ++i) {

    // Rcout << "======================================================\nProceeding to row " << i << ":\n";

    // indices of entries to be computed in diagonal band for current row
    start_ind = std::max(1, i - k); // column 0 already known
    stop_ind = std::min(m, i + k);
    // Rcout << "column indices of entries to be computed in current row: " << start_ind << " through " << stop_ind << "\n\n------------------------------------------\n\n";

    // Value above-left of leftmost entry of current band row (will be
    // sent to prev_diag before being used)
    prev_above = current_row[start_ind - 1];

    // If leftmost entry of current band row is in the leftmost column, it is
    // assigned the default value; if not, the entry to its left is treated as
    // infinite-valued
    if (i - k < 1) {
      current_row[0] = i;
    } else {
      current_row[start_ind - 1] = dist_inf;
    }

    // reset bound_exceeded
    bound_exceeded = true;

    // compute each entry of current band row
    for (int j = start_ind; j < (stop_ind + 1); ++j) {
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

    // If all entries of current band row exceed bound, we are done
    if (bound_exceeded) { return(-1); }

  }

  // If distance exceeds bound, return -1
  if (current_row[m] > k) { return(-1); }


  // return distance
  return(current_row[m]);

}