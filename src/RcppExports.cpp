// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// antibody_model_monophasic_cpp_internal
double antibody_model_monophasic_cpp_internal(int i, int t1, int b, arma::cube immune_histories, arma::cube biomarker_states, List kinetics_parameters, DataFrame biomarker_map);
RcppExport SEXP _serosim_antibody_model_monophasic_cpp_internal(SEXP iSEXP, SEXP t1SEXP, SEXP bSEXP, SEXP immune_historiesSEXP, SEXP biomarker_statesSEXP, SEXP kinetics_parametersSEXP, SEXP biomarker_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type immune_histories(immune_historiesSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type biomarker_states(biomarker_statesSEXP);
    Rcpp::traits::input_parameter< List >::type kinetics_parameters(kinetics_parametersSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type biomarker_map(biomarker_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(antibody_model_monophasic_cpp_internal(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map));
    return rcpp_result_gen;
END_RCPP
}
// antibody_model_biphasic_cpp_internal
double antibody_model_biphasic_cpp_internal(int i, int t1, int b, arma::cube immune_histories, arma::cube biomarker_states, List kinetics_parameters, DataFrame biomarker_map);
RcppExport SEXP _serosim_antibody_model_biphasic_cpp_internal(SEXP iSEXP, SEXP t1SEXP, SEXP bSEXP, SEXP immune_historiesSEXP, SEXP biomarker_statesSEXP, SEXP kinetics_parametersSEXP, SEXP biomarker_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type immune_histories(immune_historiesSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type biomarker_states(biomarker_statesSEXP);
    Rcpp::traits::input_parameter< List >::type kinetics_parameters(kinetics_parametersSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type biomarker_map(biomarker_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(antibody_model_biphasic_cpp_internal(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_serosim_antibody_model_monophasic_cpp_internal", (DL_FUNC) &_serosim_antibody_model_monophasic_cpp_internal, 7},
    {"_serosim_antibody_model_biphasic_cpp_internal", (DL_FUNC) &_serosim_antibody_model_biphasic_cpp_internal, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_serosim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
