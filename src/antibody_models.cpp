#include "utility.h"

//' Monophasic antibody boosting-waning model Rcpp implementation
//' 
//' @description Identical to \code{\link{antibody_model_monophasic}}, but implemented in Cpp
//'
//' @inheritParams antibody_model_monophasic
//'
//' @return Biomarker quantity at the specified time point
//' @export
//' @family antibody_models
//'
//' @examples
//' tmp_pars <- list()
//' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, example_model_pars_numeric)
//' antibody_model_monophasic_cpp(1,1,1,example_exposure_histories_wide, example_biomarker_states_wide, 
//' tmp_pars, example_biomarker_map_numeric)
//[[Rcpp::export]]
double antibody_model_monophasic_cpp(int i, int t1, int b,
                                     arma::cube exposure_histories,
                                     arma::cube biomarker_states,
                                     List kinetics_parameters,
                                     DataFrame biomarker_map) {
  // Find which successful exposures correspond to this biomarker
  //arma::vec exposure_id_tmp = subset_dataframe_integer(biomarker_map, "biomarker_id","exposure_id",b);
  // Get this individual's exposure history
  //arma::mat tmp_exp_history = exposure_histories.row(i-1);
  // Get only times relevant to this observation
  // Get only exposures relevant to this biomarker
  //arma::mat exp_history = get_mat_cols(tmp_exp_history.rows(0,t1-1), exposure_id_tmp-1);
  //double n_infections = sum_arma_narm(exp_history);
  
  // Set starting biomarker quantity to 0
  double biomarker_quantity = 0.0;
  DataFrame tmp_kinetics_parameters;
  if(kinetics_parameters[i-1] != R_NilValue){
    tmp_kinetics_parameters = Rcpp::as<Rcpp::DataFrame>(kinetics_parameters[i-1]);
  } else {
    return biomarker_quantity;
  }
  
  arma::vec t_infs = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","t",b, "wane");
  
  if(t_infs.n_elem > 0){
    arma::vec tmp_boost = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "boost");
    arma::vec tmp_wane = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "wane");
    
    for (int j = 0; j < tmp_boost.n_elem; j++) {
      if(t_infs[j] <= t1){
        biomarker_quantity += tmp_boost[j] * std::max(0.0, 1.0 - tmp_wane[j] * (t1 - t_infs[j]));
      }
    }
  }
  
  return biomarker_quantity;
}
//' Biphasic antibody boosting-waning model Rcpp implementation
//' 
//' @description Identical to \code{\link{antibody_model_biphasic}}, but implemented in Cpp
//'
//' @inheritParams antibody_model_monophasic
//'
//' @return Biomarker quantity at the specified time point
//' @family antibody_models
//' @export
//' @examples
//' tmp_pars <- list()
//' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, example_model_pars_numeric)
//' antibody_model_biphasic_cpp(1,1,1,example_exposure_histories_wide, example_biomarker_states_wide, 
//' tmp_pars, example_biomarker_map_numeric)
//[[Rcpp::export]]
double antibody_model_biphasic_cpp(int i, int t1, int b,
                                     arma::cube exposure_histories,
                                     arma::cube biomarker_states,
                                     List kinetics_parameters,
                                     DataFrame biomarker_map) {
  // Find which successful exposures correspond to this biomarker
  //arma::vec exposure_id_tmp = subset_dataframe_integer(biomarker_map, "biomarker_id","exposure_id",b);
  // Get this individual's exposure history
  //arma::mat tmp_exp_history = exposure_histories.row(i-1);
  // Get only times relevant to this observation
  // Get only exposures relevant to this biomarker
  //arma::mat exp_history = get_mat_cols(tmp_exp_history.rows(0,t1-1), exposure_id_tmp-1);
  //double n_infections = sum_arma_narm(exp_history);
  
  // Set starting biomarker quantity to 0
  double biomarker_quantity = 0.0;
  DataFrame tmp_kinetics_parameters;
  if(kinetics_parameters[i-1] != R_NilValue){
    tmp_kinetics_parameters = Rcpp::as<Rcpp::DataFrame>(kinetics_parameters[i-1]);
  } else {
    return biomarker_quantity;
  }
  
  arma::vec t_infs = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","t",b, "wane_long");
  
  if(t_infs.n_elem > 0){
    arma::vec tmp_boost_long = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "boost_long");
    arma::vec tmp_boost_short = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "boost_short");
    
    arma::vec tmp_wane_long = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "wane_long");
    arma::vec tmp_wane_short = subset_dataframe_numeric_twice(tmp_kinetics_parameters, "b","name","realized_value",b, "wane_short");
    
    for (int j = 0; j < t_infs.n_elem; j++) {
      if(t_infs[j] <= t1){
        biomarker_quantity += tmp_boost_long[j] * std::max(0.0, 1.0 - tmp_wane_long[j] * (t1 - t_infs[j])) + 
          tmp_boost_short[j] * std::max(0.0, 1.0 - tmp_wane_short[j] * (t1 - t_infs[j]));
      }
    }
  }
  
  return biomarker_quantity;
}

