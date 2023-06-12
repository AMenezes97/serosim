#include "utility.h"

// Internal Cpp function for monophasic antibody model
double antibody_model_monophasic_cpp_internal(int i, int t1, int b,
                                     arma::cube immune_histories,
                                     arma::cube biomarker_states,
                                     List kinetics_parameters,
                                     DataFrame biomarker_map) {
  // Find which successful exposures correspond to this biomarker
  //arma::vec exposure_id_tmp = subset_dataframe_integer(biomarker_map, "biomarker_id","exposure_id",b);
  // Get this individual's immune history
  //arma::mat tmp_imm_history = immune_histories.row(i-1);
  // Get only times relevant to this observation
  // Get only exposures relevant to this biomarker
  //arma::mat imm_history = get_mat_cols(tmp_imm_history.rows(0,t1-1), exposure_id_tmp-1);
  //double n_infections = sum_arma_narm(imm_history);
  
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

// Internal Cpp function for biphasic antibody model
double antibody_model_biphasic_cpp_internal(int i, int t1, int b,
                                     arma::cube immune_histories,
                                     arma::cube biomarker_states,
                                     List kinetics_parameters,
                                     DataFrame biomarker_map) {
  // Find which successful exposures correspond to this biomarker
  //arma::vec exposure_id_tmp = subset_dataframe_integer(biomarker_map, "biomarker_id","exposure_id",b);
  // Get this individual's immune history
  //arma::mat tmp_imm_history = immune_histories.row(i-1);
  // Get only times relevant to this observation
  // Get only exposures relevant to this biomarker
  //arma::mat imm_history = get_mat_cols(tmp_imm_history.rows(0,t1-1), exposure_id_tmp-1);
  //double n_infections = sum_arma_narm(imm_history);
  
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

