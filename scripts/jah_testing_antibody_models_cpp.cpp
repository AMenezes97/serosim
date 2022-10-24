#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

// [[Rcpp::export]]
IntegerVector cube2vector(int i, const arma::cube exposure_histories){
    IntegerVector x(exposure_histories.n_cols*exposure_histories.n_slices);
    for(int j = 0; j < x.length(); ++j){
        x[j] = exposure_histories.row(i)[j];
    }
    return x;
}

double cross_reactivity(int e, int ag, NumericVector pars){
    return 1;
}
double titer_mediation(double titer, NumericVector pars){
    return 1;
}



// [[Rcpp::export]]
double antibody_model_cpp(int i, int t, int ag, 
                          const arma::cube exposure_histories, 
                          const NumericVector pars,
                          const DataFrame antigen_map) {
    double y = 0;
    
    // 351 microseconds
    
    // Extract infection times and record exposure type
    int n_times = exposure_histories.n_cols;
    int n_exposure_types = exposure_histories.n_slices;
    
    // Get exposure types to use
    IntegerVector ags = antigen_map["Ag"];
    ags = ags - 1;
    IntegerVector es = antigen_map["E"];
    es = es - 1;
    IntegerVector es_use = es[ags == ag];
    
    // Down to 648 microseconds
    
    // Extract this individual's exposure history as a vector
    IntegerVector exp_hist = cube2vector(i, exposure_histories);
    
    // 869 microseconds
    
    IntegerVector times = rep(seq(0,exposure_histories.n_cols-1), exposure_histories.n_slices);
    IntegerVector exposure_types = rep(seq(0, exposure_histories.n_slices-1),exposure_histories.n_cols);
    
    // 1170 microseconds
    
    LogicalVector omg = in(exposure_types, es_use);
    IntegerVector exp_types = exposure_types[omg];
    
    // 1490 microseconds
    
    exp_types = exposure_types[exp_hist == 1 & times <= t & in(exposure_types,es_use)];
    NumericVector exp_times =  as<NumericVector>(times[exp_hist == 1 & times <= t & in(exposure_types,es_use)]);
    int n_exps = exp_times.length();
    
    // down to 2470 microseconds
    
    double tmp_titer = 0;
    double tmp_wane = 0;
    if(n_exps > 0){
        NumericVector boosts(n_exps);
    
        for(int j = 0; j < n_exps; ++j){
            // Get boost to this Ag prior to waning
            boosts[j] = pars["boost"] * // Unmitigated boost
                cross_reactivity(exp_types[j], ag, pars);
            
            // If this is not the first exposure
            if(j > 0){
                // Pull all boosts prior to this exposure
                NumericVector tmp_boosts = boosts[seq(0,j-1)];
                
                // Calculate what proportion of these boosts have waned by the time
                // of this exposure
                NumericVector tmp_waned = exp_times[j] - exp_times[seq(0,j-1)];
                tmp_waned = pmax(1.0 - pars["wane"]*tmp_waned,0.0);

                // Modify each boost by amount of waning
                tmp_boosts = tmp_boosts * tmp_waned;
                
                // Titer at this time is the sum of these waned boosts
                tmp_titer = sum(tmp_boosts);
                
                // Update boost amount by titer mediation
                boosts[j] = boosts[j] * titer_mediation(tmp_titer,pars);
            }
        }
        NumericVector tmp_waned = t - exp_times;
        boosts = pmax(1.0 - pars["wane"]*tmp_waned, 0.0) * boosts;

        y = sum(boosts);
    }
    return y;
}

// [[Rcpp::export]]
NumericVector antibody_model_wrapper(IntegerVector i, IntegerVector t, IntegerVector ag,
                                     const arma::cube exposure_histories, 
                                     const NumericVector pars,
                                     const DataFrame antigen_map){
    NumericVector y(i.length());
    for(int j = 0; j < i.length(); ++j){
        y[j] = antibody_model_cpp(i[j], t[j], ag[j],
                                  exposure_histories, pars, antigen_map
        );
    }
    return y;
}
