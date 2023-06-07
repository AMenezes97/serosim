#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef SUBSET_DATAFRAME_INTEGER_H
#define SUBSET_DATAFRAME_INTEGER_H
arma::vec subset_dataframe_integer(DataFrame mydata, String condition_vector, String return_vector, int condition);
#endif

#ifndef SUBSET_DATAFRAME_NUMERIC_TWICE_H
#define SUBSET_DATAFRAME_NUMERIC_TWICE_H
arma::vec subset_dataframe_numeric_twice(DataFrame x, 
                                         String condition_vector1, 
                                         String condition_vector2,
                                         String return_vector, 
                                         double condition1, 
                                         String condition2);
#endif


#ifndef SUBSET_DATAFRAME_NUMERIC_H
#define SUBSET_DATAFRAME_NUMERIC_H
arma::vec subset_dataframe_numeric(DataFrame mydata, String condition_vector, String return_vector, double condition);
#endif

#ifndef SUBSET_DATAFRAME_CHAR_H
#define SUBSET_DATAFRAME_CHAR_H
arma::vec subset_dataframe_char(DataFrame mydata, String condition_vector, String return_vector, String condition);
#endif

#ifndef GET_CUBE_SLICES_H
#define GET_CUBE_SLICES_H
arma::cube get_cube_slices(arma::cube x, arma::vec indices);
#endif

#ifndef GET_MAT_COLS_H
#define GET_MAT_COLS_H
arma::mat get_mat_cols(arma::mat x, arma::vec indices);
#endif

#ifndef SUM_ARMA_NA_RM_H
#define SUM_ARMA_NA_RM_H
double sum_arma_na_rm(arma::mat& X);
#endif
