#' Biphasic antibody boosting-waning model
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @param i Individual
#' @param t1 time
#' @param ag antigen
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states An array of antibody states across all individuals, time steps and antigen IDs
#' @param kinetics_parameters An object of all kinetics parameters for all exposures
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param ... 
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_biphasic <-  function(i, t1, ag, exposure_histories, antibody_states, kinetics_parameters, antigen_map, ...){
  ## Find which successful exposures correspond to this antigen 
  exposure_id_tmp<-antigen_map$exposure_id[antigen_map$antigen_id==ag]
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting titer to 0
  titer<-0
  
  ## Calculate current titer if there has been an exposure 
  if(sum(exp_history,na.rm = TRUE)==0){
    return(0)
  }
  if(sum(exp_history,na.rm = TRUE)>0){
    
    ## Extract all kinetics_parameters for antigen 
    ag_tmp<-ag
    
    tmp_kinetics_parameters <- data.table(kinetics_parameters[[i]])
    tmp_kinetics_parameters<-tmp_kinetics_parameters[t<t1 & ag==ag_tmp]
    
    tmp_boost_long <- tmp_kinetics_parameters[name == "boost_long"] 
    tmp_boost_short <- tmp_kinetics_parameters[name == "boost_short"] 
    
    tmp_wane_long <- tmp_kinetics_parameters[name == "wane_long"] 
    tmp_wane_short <- tmp_kinetics_parameters[name == "wane_short"] 
    
    
    for(j in seq_along(tmp_boost_long$realized_value)){
      titer<- titer + tmp_boost_long$realized_value[j]*max(0,1-tmp_wane_long$realized_value[j]*(t1-tmp_wane_long$t[j])) 
      + tmp_boost_short$realized_value[j]*max(0,1-tmp_wane_short$realized_value[j]*(t1-tmp_wane_short$t[j]))
    }
    titer
    
  }
}

