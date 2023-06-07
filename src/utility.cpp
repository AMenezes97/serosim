
#include "utility.h"

// [[Rcpp::export]]
arma::vec subset_dataframe_integer(DataFrame mydata, String condition_vector, String return_vector, int condition) {
  IntegerVector vector1 = mydata[condition_vector];
  NumericVector vector2 = mydata[return_vector];
  NumericVector result;
  
  for (int i = 0; i < vector1.size(); i++) {
    if (vector1[i] == condition) {
      result.push_back(vector2[i]);
    }
  }
  return as<arma::vec>(wrap(result));
}

// [[Rcpp::export]]
arma::vec subset_dataframe_numeric_twice(DataFrame x, 
                                             String condition_vector1, 
                                             String condition_vector2,
                                             String return_vector, 
                                             double condition1, 
                                             String condition2){
  IntegerVector vector1 = x[condition_vector1];
  CharacterVector vector2 = x[condition_vector2];
  NumericVector vector3 = x[return_vector];
  NumericVector result;
  
  for (int i = 0; i < vector1.size(); i++) {
    if (vector1[i] == condition1 && vector2[i] == condition2) {
      result.push_back(vector3[i]);
    }
  }
  return as<arma::vec>(wrap(result));
}

// [[Rcpp::export]]
arma::vec subset_dataframe_numeric(DataFrame mydata, String condition_vector, String return_vector, double condition) {
  NumericVector vector1 = mydata[condition_vector];
  NumericVector vector2 = mydata[return_vector];
  
  NumericVector result;
  
  for (int i = 0; i < vector1.size(); i++) {
    if (vector1[i] == condition) {
      result.push_back(vector2[i]);
    }
  }
  return as<arma::vec>(wrap(result));
}

// [[Rcpp::export]]
arma::vec subset_dataframe_char(DataFrame mydata, String condition_vector, String return_vector, String condition) {
  CharacterVector vector1 = mydata[condition_vector];
  NumericVector vector2 = mydata[return_vector];
  
  NumericVector result;
  
  for (int i = 0; i < vector1.size(); i++) {
    if (vector1[i] == condition) {
      result.push_back(vector2[i]);
    }
  }
  return as<arma::vec>(wrap(result));
}
// [[Rcpp::export]]
arma::cube get_cube_slices(arma::cube x, arma::vec indices){
  arma::cube y(x.n_rows,x.n_cols, indices.n_elem);
  for(int i=0; i<indices.n_elem; ++i){
    y.slice(i) = x.slice(indices(i));
  }
  return y;
}
// [[Rcpp::export]]
arma::mat get_mat_cols(arma::mat x, arma::vec indices){
  arma::mat y(x.n_rows,indices.n_elem);
  for(int i=0; i<indices.n_elem; ++i){
    y.col(i) = x.col(indices(i));
  }
  return y;
}


// [[Rcpp::export]]
double sum_arma_na_rm(arma::mat& X) {
  double sum = 0;
  for (int i = 0; i < X.size(); ++i) {
    if (arma::is_finite(X(i)))
      sum += X(i);
  }
  return sum;
}

