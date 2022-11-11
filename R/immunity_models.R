#' Immunity Model Where All Exposures Are Successful
#'
#' @description Probability of success is 1 so all exposure events are successful
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_all_successful <- function(i, t, x, exposure_histories, 
                           antibody_states, demography, biomarker_map, model_pars, ...){
  return(1)
}

#' Immunity Model For Vaccination Events Only
#' 
#' @description This immunity model should only be used if all exposures are vaccination events. The probability of successful exposure(vaccination event) depends on the number of vaccines an individual has received prior to time t. If the individual is under the maximum vaccinations allotted then the probability of successful exposure is 1.
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_vacc_events A vector of the maximum number of vaccination events possible for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_only <- function(i, t, x, exposure_histories, 
                              antibody_states, demography, biomarker_map, 
                              model_pars, max_vacc_events, vacc_age,...){
  ## Calculate the individual's current age
  birth_time<-unique(demography$birth[demography$i==i])
  curr_age<- t-birth_time
  
  ## If the individual is above the minimum age of vaccination 
  if(curr_age>=vacc_age[x]){
  ## Count the total number of successful exposures to e thus far 
    curr_vacc_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
    ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
    if(curr_vacc_events<max_vacc_events[x]){
      return(1)
    }else{
      return(0)
    }
  }else{
    return(0)
  }
}
 

#' Simple Immunity Model For Vaccination Events and Natural Infection Events
#' 
#' @description This immunity model should only be used with vaccines and natural infection. The probability of successful exposure for vaccination events depends on the number of vaccines an individual has received prior to time t and their current age. If the individual is under the maximum vaccinations allotted and is of an age eligible for vaccination then the probability of successful exposure is 1. The probability of a successful natural infection is dependent on the total number of infections that individual has experienced thus far. If the individual is under the maximum number of infections allotted then the probability of successful exposure is 1. 
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_events A vector of the maximum number of successful exposure events possible for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param vacc_exposures A vector of exposure IDs (x) which represents vaccination events
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_ifxn_simple <- function(i, t, x, exposure_histories, 
                                     antibody_states, demography, biomarker_map, model_pars, max_events, vacc_exposures, vacc_age, ...){
  ## If an exposure event is a vaccination event, then guaranteed exposure unless the individual has already been vaccinated
  if(x %in% c(vacc_exposures)){  	
    ## Calculate the individual's current age
    birth_time<-unique(demography$birth[demography$i==i])
    curr_age<- t-birth_time
    
    ## If the individual is above the minimum age of vaccination 
    if(curr_age>=vacc_age[x]){
         ## Count the total number of successful exposures to e thus far 
          curr_vacc_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
        ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
        if(curr_vacc_events<max_events[x]){
           return(1)
             }else{
             return(0)
              }
    }
    if(curr_age<vacc_age[x]){
      return(0)
    }
  }
  else{ ## If the exposure event is an infection event
    ## Count the total number of successful exposures to e thus far 
    curr_ifx_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
    if(curr_ifx_events<max_events[x]){
      return(1)
    }else{
      return(0)
  }
  }
}



#' Immunity Model For Natural Infection Events With Titer-Mediated Protection
#' 
#' @description  This immunity model should only be used if all exposures are natural infection events. The probability of successful exposure is dependent on the individual’s antibody titer at the time of exposure. 
#'    
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_ifxn_titer_prot <- function(i, t, x, exposure_histories, 
                              antibody_states, demography, biomarker_map, model_pars, ...){
    ## Find biomarkers which are boosted by this exposure type
    ## The assumption here is that the titer levels to these biomarkers will determine if an individual is protected
    b<-biomarker_map$biomarker_id[biomarker_map$exposure_id==x]
    ## Find current titer to all relevant biomarkers
    curr_t <- antibody_states[i,t,b] ## How to deal with titers against multiple biomarkers? Should they be added?
    curr_t <- sum(curr_t)
    
    ## Pull out necessary variables 
    titer_prot_midpoint <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="titer_prot_midpoint"] ## Would each biomarker have it's own value?
    titer_prot_width <- model_pars$mean[model_pars$exposure_id==x  & model_pars$name=="titer_prot_width"]
    
    ## Create a function to calculate the risk of infection at a given titer
    titer_protection <- function(titer, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(titer - alpha1)))
      return(risk)
    }
    
    prob_success<- (1-titer_protection(curr_t, titer_prot_midpoint, titer_prot_width))
    
    return(prob_success)
  }


#' Immunity Model For Vaccination Events and Natural Infection Events With Titer-Mediated Protection
#' 
#' @description This immunity model should be used if exposures represent vaccination and natural infection events. The probability of successful vaccination exposure depends on the number of vaccines received prior to time t while the probability of successful infection is dependent on the titer at the time of exposure.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_vacc_events A vector of the maximum number of vaccination events possible for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param vacc_exposures A vector of exposure IDs (x) which represents vaccination events
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_ifxn_titer_prot <- function(i, t, x, exposure_histories,
                           antibody_states, demography, biomarker_map, model_pars, max_vacc_events, vacc_exposures, vacc_age=1, ...){
  ## If an exposure event is a vaccination event, then guaranteed exposure unless the individual has already been vaccinated
  if(x %in% c(vacc_exposures)){  
    ## Calculate the individual's current age
    birth_time<-unique(demography$birth[demography$i==i])
    curr_age<- t-birth_time
    
    ## If the individual is above the minimum age of vaccination 
    if(curr_age>=vacc_age[x]){
    ## Count the total number of successful exposures to e thus far 
    curr_vacc_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
    ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
    if(curr_vacc_events<max_vacc_events[x]){
      return(1)
    }else{
      return(0)
    }
    }
    if(curr_age<vacc_age[x]){
      return(0)
    }
  } 
  else {
    ## Find biomarkers which are boosted by this exposure type
    ## The assumption here is that the titer levels to these biomarkers will determine if an individual is protected
    b<-biomarker_map$biomarker_id[biomarker_map$exposure_id==x]
    ## Find current titer to all relevant biomarkers
    curr_t <- antibody_states[i,t,b] ## How to deal with titers against multiple biomarkers? Should they be added?
    curr_t <- sum(curr_t)
    
    ## Pull out necessary variables 
    titer_prot_midpoint <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="titer_prot_midpoint"]
    titer_prot_width <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="titer_prot_width"]
    
    ## Create a function to calculate the risk of infection at a given titer
    titer_protection <- function(titer, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(titer - alpha1)))
      return(risk)
    }
    
    prob_success<- (1-titer_protection(curr_t, titer_prot_midpoint, titer_prot_width))
    
    return(prob_success)
  }
}

