// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hamAdjacencyMatSparse
arma::sp_umat hamAdjacencyMatSparse(std::vector<std::string> strings, const int& maxdist);
RcppExport SEXP _RepSeqNetworkAnalysis_hamAdjacencyMatSparse(SEXP stringsSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type strings(stringsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxdist(maxdistSEXP);
    rcpp_result_gen = Rcpp::wrap(hamAdjacencyMatSparse(strings, maxdist));
    return rcpp_result_gen;
END_RCPP
}
// hamDistBounded
int hamDistBounded(std::string a, std::string b, const int& k);
RcppExport SEXP _RepSeqNetworkAnalysis_hamDistBounded(SEXP aSEXP, SEXP bSEXP, SEXP kSEXP) {
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
arma::sp_umat levAdjacencyMatSparse(std::vector<std::string> strings, const int& maxdist);
RcppExport SEXP _RepSeqNetworkAnalysis_levAdjacencyMatSparse(SEXP stringsSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type strings(stringsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxdist(maxdistSEXP);
    rcpp_result_gen = Rcpp::wrap(levAdjacencyMatSparse(strings, maxdist));
    return rcpp_result_gen;
END_RCPP
}
// levDistBounded
int levDistBounded(std::string a, std::string b, const int& k);
RcppExport SEXP _RepSeqNetworkAnalysis_levDistBounded(SEXP aSEXP, SEXP bSEXP, SEXP kSEXP) {
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
    {"_RepSeqNetworkAnalysis_hamAdjacencyMatSparse", (DL_FUNC) &_RepSeqNetworkAnalysis_hamAdjacencyMatSparse, 2},
    {"_RepSeqNetworkAnalysis_hamDistBounded", (DL_FUNC) &_RepSeqNetworkAnalysis_hamDistBounded, 3},
    {"_RepSeqNetworkAnalysis_levAdjacencyMatSparse", (DL_FUNC) &_RepSeqNetworkAnalysis_levAdjacencyMatSparse, 2},
    {"_RepSeqNetworkAnalysis_levDistBounded", (DL_FUNC) &_RepSeqNetworkAnalysis_levDistBounded, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RepSeqNetworkAnalysis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}