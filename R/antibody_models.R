## What's t1 again? t+1?
## Since the antibody model is called within Ag loop; the model should be structured to look at one antigen at a time? 

#' Antibody Model Version 1
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @param i Individual
#' @param t time
#' @param ag antigen
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param kinetics_parameters An object of all kinetics parameters for all exposures
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_V1 <- function(i, t, ag, exposure_histories, kinetics_parameters, antigen_map){
  ## Find which successful exposures correspond to this antigen 
  antigen_map_tmp<- antigen_map %>%  filter(antigen_id==ag)
  exposure_id_tmp<-antigen_map_tmp$exposure_id
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:(t1-1),exposure_id_tmp]
  
  ## Set starting titer to 0
  titer<-0
  
  t1<- t+1 ## CHeck this with James! 
  
  ## Calculate current titer if there has been an exposure 
  if(t1 > 1 & sum(exp_history > 0)){
    
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
    

  }
 ## Returns current titer
  titer
}
