// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ordinal_bart_main
Rcpp::List ordinal_bart_main(arma::ivec const& y, arma::mat const& x, arma::mat const& xtest, int burn, int nd, int m, int nc, double pbd, double pb, double alpha, double betap, double kappa, arma::mat const& Ad, double s, arma::mat const& inc_root, arma::vec const& dstarbar);
RcppExport SEXP _OrderedBart_ordinal_bart_main(SEXP ySEXP, SEXP xSEXP, SEXP xtestSEXP, SEXP burnSEXP, SEXP ndSEXP, SEXP mSEXP, SEXP ncSEXP, SEXP pbdSEXP, SEXP pbSEXP, SEXP alphaSEXP, SEXP betapSEXP, SEXP kappaSEXP, SEXP AdSEXP, SEXP sSEXP, SEXP inc_rootSEXP, SEXP dstarbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type xtest(xtestSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< double >::type pbd(pbdSEXP);
    Rcpp::traits::input_parameter< double >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type betap(betapSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Ad(AdSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type inc_root(inc_rootSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type dstarbar(dstarbarSEXP);
    rcpp_result_gen = Rcpp::wrap(ordinal_bart_main(y, x, xtest, burn, nd, m, nc, pbd, pb, alpha, betap, kappa, Ad, s, inc_root, dstarbar));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OrderedBart_ordinal_bart_main", (DL_FUNC) &_OrderedBart_ordinal_bart_main, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_OrderedBart(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}