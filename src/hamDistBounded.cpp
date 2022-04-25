//compute Hamming distance between two strings a and b of lengths n and m,
//subject to a specified upper bound k.  If a and b differ in length, the longer
//string is truncated to the length of the shorter string and the difference in
//length is added to the distance.  If it is found that the Hamming distance
//will exceed the bound, the computation is halted and a value of -1 is
//returned.

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <strings.h>
using namespace Rcpp;

// [[Rcpp::export]]
int hamDistBounded(std::string a,
                   std::string b,
                   const int& k) {

  // if strings are identical, return 0
  if (a == b) { return(0); }

  // get string lengths
  int n = a.length();
  int m = b.length();
  // Rcout << "length of a: " << n << " characters\n";
  // Rcout << "length of b: " << m << " characters\n";

  // Initialize distance value, bound below by difference in string length
  int dist = abs(n - m);

  // stop if distance exceeds bound
  if (dist > k) { return(-1); }

  // Compute hamming distance; longer string truncated to length of shorter
  for (int i = 0; i < std::min(n, m); ++i) {
    if (a[i] != b[i]) { dist++; }
    if (dist > k) { return(-1); } // stop if distance exceeds bound
  }

  // return distance
  return(dist);

}