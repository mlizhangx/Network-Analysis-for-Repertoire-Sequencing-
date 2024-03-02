#ifndef DROPDEGREEZERO_H
#define DROPDEGREEZERO_H

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

void dropDegreeZero(
    bool drop_deg_zero,
    arma::sp_umat& out,
    std::string tempfile
);

#endif  /* DROPDEGREEZERO_H */