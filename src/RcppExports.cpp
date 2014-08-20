// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mmvnpdfC
NumericMatrix mmvnpdfC(NumericMatrix x, NumericMatrix mean, List varcovM);
RcppExport SEXP NPflow_mmvnpdfC(SEXP xSEXP, SEXP meanSEXP, SEXP varcovMSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP );
        Rcpp::traits::input_parameter< List >::type varcovM(varcovMSEXP );
        NumericMatrix __result = mmvnpdfC(x, mean, varcovM);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mmvsnpdfC
NumericMatrix mmvsnpdfC(NumericMatrix x, NumericMatrix xi, NumericMatrix psi, List sigma);
RcppExport SEXP NPflow_mmvsnpdfC(SEXP xSEXP, SEXP xiSEXP, SEXP psiSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type xi(xiSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type psi(psiSEXP );
        Rcpp::traits::input_parameter< List >::type sigma(sigmaSEXP );
        NumericMatrix __result = mmvsnpdfC(x, xi, psi, sigma);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mvnpdfC
NumericVector mvnpdfC(NumericMatrix x, NumericVector mean, NumericMatrix varcovM);
RcppExport SEXP NPflow_mvnpdfC(SEXP xSEXP, SEXP meanSEXP, SEXP varcovMSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type varcovM(varcovMSEXP );
        NumericVector __result = mvnpdfC(x, mean, varcovM);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
