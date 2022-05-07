// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// permuC
arma::vec permuC(arma::vec x);
RcppExport SEXP _LA_permuC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(permuC(x));
    return rcpp_result_gen;
END_RCPP
}
// seqC
arma::vec seqC(unsigned a, int b);
RcppExport SEXP _LA_seqC(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(seqC(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rLHDC
arma::mat rLHDC(int n, int k);
RcppExport SEXP _LA_rLHDC(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(rLHDC(n, k));
    return rcpp_result_gen;
END_RCPP
}
// dijC
double dijC(arma::mat X, int i, int j, int q);
RcppExport SEXP _LA_dijC(SEXP XSEXP, SEXP iSEXP, SEXP jSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(dijC(X, i, j, q));
    return rcpp_result_gen;
END_RCPP
}
// phi_pC
double phi_pC(arma::mat X, int p, int q);
RcppExport SEXP _LA_phi_pC(SEXP XSEXP, SEXP pSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_pC(X, p, q));
    return rcpp_result_gen;
END_RCPP
}
// MaxProCriterionC
double MaxProCriterionC(arma::mat X);
RcppExport SEXP _LA_MaxProCriterionC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MaxProCriterionC(X));
    return rcpp_result_gen;
END_RCPP
}
// corC
double corC(arma::vec x, arma::vec y);
RcppExport SEXP _LA_corC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(corC(x, y));
    return rcpp_result_gen;
END_RCPP
}
// MaxAbsCorC
double MaxAbsCorC(arma::mat X);
RcppExport SEXP _LA_MaxAbsCorC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MaxAbsCorC(X));
    return rcpp_result_gen;
END_RCPP
}
// AvgAbsCorC
double AvgAbsCorC(arma::mat X);
RcppExport SEXP _LA_AvgAbsCorC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(AvgAbsCorC(X));
    return rcpp_result_gen;
END_RCPP
}
// exchangeC
arma::mat exchangeC(arma::mat X, int j, String type);
RcppExport SEXP _LA_exchangeC(SEXP XSEXP, SEXP jSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(exchangeC(X, j, type));
    return rcpp_result_gen;
END_RCPP
}
// LA_LHDC
arma::mat LA_LHDC(int n, int k, int m, int N, String OC, int p, int q);
RcppExport SEXP _LA_LA_LHDC(SEXP nSEXP, SEXP kSEXP, SEXP mSEXP, SEXP NSEXP, SEXP OCSEXP, SEXP pSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< String >::type OC(OCSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(LA_LHDC(n, k, m, N, OC, p, q));
    return rcpp_result_gen;
END_RCPP
}
// factorialC
double factorialC(int x);
RcppExport SEXP _LA_factorialC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(factorialC(x));
    return rcpp_result_gen;
END_RCPP
}
// modC
int modC(int a, int b);
RcppExport SEXP _LA_modC(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(modC(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rOofAC
arma::mat rOofAC(int n, int k);
RcppExport SEXP _LA_rOofAC(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(rOofAC(n, k));
    return rcpp_result_gen;
END_RCPP
}
// PWOC
arma::mat PWOC(arma::mat X);
RcppExport SEXP _LA_PWOC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(PWOC(X));
    return rcpp_result_gen;
END_RCPP
}
// TC
arma::mat TC(arma::mat X);
RcppExport SEXP _LA_TC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(TC(X));
    return rcpp_result_gen;
END_RCPP
}
// MOMC
arma::mat MOMC(arma::mat X);
RcppExport SEXP _LA_MOMC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MOMC(X));
    return rcpp_result_gen;
END_RCPP
}
// LA_OofAC
arma::mat LA_OofAC(int n, int k, int m, int N);
RcppExport SEXP _LA_LA_OofAC(SEXP nSEXP, SEXP kSEXP, SEXP mSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(LA_OofAC(n, k, m, N));
    return rcpp_result_gen;
END_RCPP
}
// D
double D(arma::mat X);
RcppExport SEXP _LA_D(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(D(X));
    return rcpp_result_gen;
END_RCPP
}
// A
double A(arma::mat X);
RcppExport SEXP _LA_A(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(A(X));
    return rcpp_result_gen;
END_RCPP
}
// GscoreC
double GscoreC(arma::mat X, arma::vec x);
RcppExport SEXP _LA_GscoreC(SEXP XSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(GscoreC(X, x));
    return rcpp_result_gen;
END_RCPP
}
// rSign
double rSign(int m);
RcppExport SEXP _LA_rSign(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rSign(m));
    return rcpp_result_gen;
END_RCPP
}
// G
double G(arma::mat Y);
RcppExport SEXP _LA_G(SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(G(Y));
    return rcpp_result_gen;
END_RCPP
}
// LA_OptC
arma::mat LA_OptC(int n, arma::vec lb, arma::vec ub, int m, int N, String OC, double alpha);
RcppExport SEXP _LA_LA_OptC(SEXP nSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP mSEXP, SEXP NSEXP, SEXP OCSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< String >::type OC(OCSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(LA_OptC(n, lb, ub, m, N, OC, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LA_permuC", (DL_FUNC) &_LA_permuC, 1},
    {"_LA_seqC", (DL_FUNC) &_LA_seqC, 2},
    {"_LA_rLHDC", (DL_FUNC) &_LA_rLHDC, 2},
    {"_LA_dijC", (DL_FUNC) &_LA_dijC, 4},
    {"_LA_phi_pC", (DL_FUNC) &_LA_phi_pC, 3},
    {"_LA_MaxProCriterionC", (DL_FUNC) &_LA_MaxProCriterionC, 1},
    {"_LA_corC", (DL_FUNC) &_LA_corC, 2},
    {"_LA_MaxAbsCorC", (DL_FUNC) &_LA_MaxAbsCorC, 1},
    {"_LA_AvgAbsCorC", (DL_FUNC) &_LA_AvgAbsCorC, 1},
    {"_LA_exchangeC", (DL_FUNC) &_LA_exchangeC, 3},
    {"_LA_LA_LHDC", (DL_FUNC) &_LA_LA_LHDC, 7},
    {"_LA_factorialC", (DL_FUNC) &_LA_factorialC, 1},
    {"_LA_modC", (DL_FUNC) &_LA_modC, 2},
    {"_LA_rOofAC", (DL_FUNC) &_LA_rOofAC, 2},
    {"_LA_PWOC", (DL_FUNC) &_LA_PWOC, 1},
    {"_LA_TC", (DL_FUNC) &_LA_TC, 1},
    {"_LA_MOMC", (DL_FUNC) &_LA_MOMC, 1},
    {"_LA_LA_OofAC", (DL_FUNC) &_LA_LA_OofAC, 4},
    {"_LA_D", (DL_FUNC) &_LA_D, 1},
    {"_LA_A", (DL_FUNC) &_LA_A, 1},
    {"_LA_GscoreC", (DL_FUNC) &_LA_GscoreC, 2},
    {"_LA_rSign", (DL_FUNC) &_LA_rSign, 1},
    {"_LA_G", (DL_FUNC) &_LA_G, 1},
    {"_LA_LA_OptC", (DL_FUNC) &_LA_LA_OptC, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_LA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
