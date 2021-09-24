// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// grad_fun_cpp
NumericVector grad_fun_cpp(NumericVector nlm_par, DataFrame dset, double lambda_theta, double lambda_in, double lambda_out, double lambda_delta, CharacterVector est_par, CharacterVector fix_par, IntegerVector est_length, IntegerVector fix_length, NumericVector fix_value, NumericMatrix groups_map, IntegerVector mt_vek, IntegerVector mt_idx, IntegerVector dim_resp, IntegerVector all_cat, IntegerVector n_th, IntegerVector XN, IntegerVector XNA, double eps, bool is_penalized_gamma, bool is_penalized_theta, bool is_penalized_delta);
RcppExport SEXP _autoRasch_grad_fun_cpp(SEXP nlm_parSEXP, SEXP dsetSEXP, SEXP lambda_thetaSEXP, SEXP lambda_inSEXP, SEXP lambda_outSEXP, SEXP lambda_deltaSEXP, SEXP est_parSEXP, SEXP fix_parSEXP, SEXP est_lengthSEXP, SEXP fix_lengthSEXP, SEXP fix_valueSEXP, SEXP groups_mapSEXP, SEXP mt_vekSEXP, SEXP mt_idxSEXP, SEXP dim_respSEXP, SEXP all_catSEXP, SEXP n_thSEXP, SEXP XNSEXP, SEXP XNASEXP, SEXP epsSEXP, SEXP is_penalized_gammaSEXP, SEXP is_penalized_thetaSEXP, SEXP is_penalized_deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nlm_par(nlm_parSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type dset(dsetSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_theta(lambda_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_in(lambda_inSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_out(lambda_outSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_delta(lambda_deltaSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type est_par(est_parSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type fix_par(fix_parSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type est_length(est_lengthSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type fix_length(fix_lengthSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fix_value(fix_valueSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type groups_map(groups_mapSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mt_vek(mt_vekSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mt_idx(mt_idxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dim_resp(dim_respSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type all_cat(all_catSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_th(n_thSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XN(XNSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XNA(XNASEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type is_penalized_gamma(is_penalized_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type is_penalized_theta(is_penalized_thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type is_penalized_delta(is_penalized_deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_fun_cpp(nlm_par, dset, lambda_theta, lambda_in, lambda_out, lambda_delta, est_par, fix_par, est_length, fix_length, fix_value, groups_map, mt_vek, mt_idx, dim_resp, all_cat, n_th, XN, XNA, eps, is_penalized_gamma, is_penalized_theta, is_penalized_delta));
    return rcpp_result_gen;
END_RCPP
}
// ll_cpp
double ll_cpp(arma::vec theta, arma::vec gamma, arma::mat delta, arma::mat groups, arma::mat beta, arma::vec m_cat, arma::mat X, bool gamma_penalized, bool delta_penalized, bool theta_penalized, double lambda_in, double lambda_out, double lambda_delta, double lambda_theta, double eps);
RcppExport SEXP _autoRasch_ll_cpp(SEXP thetaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP m_catSEXP, SEXP XSEXP, SEXP gamma_penalizedSEXP, SEXP delta_penalizedSEXP, SEXP theta_penalizedSEXP, SEXP lambda_inSEXP, SEXP lambda_outSEXP, SEXP lambda_deltaSEXP, SEXP lambda_thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m_cat(m_catSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type gamma_penalized(gamma_penalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type delta_penalized(delta_penalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type theta_penalized(theta_penalizedSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_in(lambda_inSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_out(lambda_outSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_delta(lambda_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_theta(lambda_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ll_cpp(theta, gamma, delta, groups, beta, m_cat, X, gamma_penalized, delta_penalized, theta_penalized, lambda_in, lambda_out, lambda_delta, lambda_theta, eps));
    return rcpp_result_gen;
END_RCPP
}
// grad_cpp
Rcpp::List grad_cpp(arma::vec theta, arma::vec gamma, arma::mat delta, arma::mat groups, arma::mat beta, arma::vec m_cat, arma::mat X, bool gamma_penalized, bool delta_penalized, bool theta_penalized, double lambda_in, double lambda_out, double lambda_delta, double lambda_theta, double eps);
RcppExport SEXP _autoRasch_grad_cpp(SEXP thetaSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP groupsSEXP, SEXP betaSEXP, SEXP m_catSEXP, SEXP XSEXP, SEXP gamma_penalizedSEXP, SEXP delta_penalizedSEXP, SEXP theta_penalizedSEXP, SEXP lambda_inSEXP, SEXP lambda_outSEXP, SEXP lambda_deltaSEXP, SEXP lambda_thetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m_cat(m_catSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type gamma_penalized(gamma_penalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type delta_penalized(delta_penalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type theta_penalized(theta_penalizedSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_in(lambda_inSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_out(lambda_outSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_delta(lambda_deltaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_theta(lambda_thetaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_cpp(theta, gamma, delta, groups, beta, m_cat, X, gamma_penalized, delta_penalized, theta_penalized, lambda_in, lambda_out, lambda_delta, lambda_theta, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_autoRasch_grad_fun_cpp", (DL_FUNC) &_autoRasch_grad_fun_cpp, 23},
    {"_autoRasch_ll_cpp", (DL_FUNC) &_autoRasch_ll_cpp, 15},
    {"_autoRasch_grad_cpp", (DL_FUNC) &_autoRasch_grad_cpp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_autoRasch(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
