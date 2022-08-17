#' Biphasic antibody boosting-waning model
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @param i Individual
#' @param t1 time
#' @param ag antigen
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param kinetics_parameters An object of all kinetics parameters for all exposures
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param antibody_states 
#' @param ... 
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_biphasic <- function(i, t1, ag, exposure_histories, antibody_states, kinetics_parameters, antigen_map, ...){
  ## Find which successful exposures correspond to this antigen 
  antigen_map_tmp<- antigen_map %>%  filter(antigen_id==ag)
  exposure_id_tmp<-antigen_map_tmp$exposure_id
  
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
    
    tmp_boost_long <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "boost_long" & ag == ag_tmp) 
    tmp_boost_short <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "boost_short" & ag == ag_tmp) 
      
    tmp_wane_long <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "wane_long" & ag == ag_tmp) 
    tmp_wane_short <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "wane_short" & ag == ag_tmp) 
    
    
    for(j in seq_along(tmp_boost_long$value)){
      titer<- titer + tmp_boost_long$value[j]*max(0,1-tmp_wane_long$value[j]*(t1-tmp_wane_long$t[j])) 
              + tmp_boost_short$value[j]*max(0,1-tmp_wane_short$value[j]*(t1-tmp_wane_short$t[j]))
    }
    titer

  }
}

#' Biphasic antibody boosting-waning model with titer-dependent boosting
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @param i Individual
#' @param t1 
#' @param ag antigen
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param kinetics_parameters An object of all kinetics parameters for all exposures
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param ... 
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_biphasic_titer_dep <- function(i, t1, ag, exposure_histories, kinetics_parameters, antigen_map, ...){
  ## Find which successful exposures correspond to this antigen 
  antigen_map_tmp<- antigen_map %>%  filter(antigen_id==ag)
  exposure_id_tmp<-antigen_map_tmp$exposure_id
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting titer to 0
  titer<-0
  
  ## Calculate current titer if there has been an exposure 
  if(is.na(sum(exp_history))){
    return(0)
  }
  if(t1 > 1 & sum(exp_history > 0)){
    
    ## Extract all kinetics_parameters for antigen 
    ag_tmp<-ag
    
    tmp_boost_long <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "boost_long" & ag == ag_tmp) 
    tmp_boost_short <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "boost_short" & ag == ag_tmp) 
    
    tmp_wane_long <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "wane_long" & ag == ag_tmp) 
    tmp_wane_short <- kinetics_parameters[[i]] %>% filter(t < t1) %>% filter(name == "wane_short" & ag == ag_tmp) 
    
    ## Goal: Adjust all boost_long and boost_short to take into account titer at time of boosting 
    ## This becomes complicated because we are currently working in the Ag level but we 
    ## might need to consider other Ag when adjusting the boost levels.
    ## For ex: The boosts in tmp_boost_long might be coming from different types of exposures 
    ## were titer levels to multiple Ags might affect titer ceiling effects.
    
    ## What if we kept titer-dependent boosting in draw_parameters(where is originally was)
    ## but we saved the true boost and realized boost.
    ## This way kinetics_parameters[[i]] will have both values so we can pull the realized 
    ## boost for the antibody_kinetics model. 
    
    
    for(j in seq_along(tmp_boost_long$value)){
      titer<- titer + tmp_boost_long$value[j]*max(0,1-tmp_wane_long$value[j]*(t1-tmp_wane_long$t[j])) 
      + tmp_boost_short$value[j]*max(0,1-tmp_wane_short$value[j]*(t1-tmp_wane_short$t[j]))
    }
    ## Returns current titer
    titer
    
  }
  else{
    return(0)
  }
}
