
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hamAdjacencyMatSparse
arma::sp_umat hamAdjacencyMatSparse(std::vector<std::string> strings, const int& maxdist, bool drop_deg_zero);
RcppExport SEXP _NAIR_hamAdjacencyMatSparse(SEXP stringsSEXP, SEXP maxdistSEXP, SEXP drop_deg_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type strings(stringsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxdist(maxdistSEXP);
    Rcpp::traits::input_parameter< bool >::type drop_deg_zero(drop_deg_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(hamAdjacencyMatSparse(strings, maxdist, drop_deg_zero));
    return rcpp_result_gen;
END_RCPP
}
// hamDistBounded
int hamDistBounded(std::string a, std::string b, const int& k);
RcppExport SEXP _NAIR_hamDistBounded(SEXP aSEXP, SEXP bSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(hamDistBounded(a, b, k));
    return rcpp_result_gen;
END_RCPP
}
// levAdjacencyMatSparse
arma::sp_umat levAdjacencyMatSparse(std::vector<std::string> strings, const int& maxdist, bool drop_deg_zero);
RcppExport SEXP _NAIR_levAdjacencyMatSparse(SEXP stringsSEXP, SEXP maxdistSEXP, SEXP drop_deg_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type strings(stringsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxdist(maxdistSEXP);
    Rcpp::traits::input_parameter< bool >::type drop_deg_zero(drop_deg_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(levAdjacencyMatSparse(strings, maxdist, drop_deg_zero));
    return rcpp_result_gen;
END_RCPP
}
// levDistBounded
int levDistBounded(std::string a, std::string b, const int& k);
RcppExport SEXP _NAIR_levDistBounded(SEXP aSEXP, SEXP bSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(levDistBounded(a, b, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NAIR_hamAdjacencyMatSparse", (DL_FUNC) &_NAIR_hamAdjacencyMatSparse, 3},
    {"_NAIR_hamDistBounded", (DL_FUNC) &_NAIR_hamDistBounded, 3},
    {"_NAIR_levAdjacencyMatSparse", (DL_FUNC) &_NAIR_levAdjacencyMatSparse, 3},
    {"_NAIR_levDistBounded", (DL_FUNC) &_NAIR_levDistBounded, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_NAIR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
