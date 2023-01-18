#' Immunity Model Where All Exposures Are Successful
#'
#' @description Probability of success is 1 so all exposure events are successful
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure events
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs
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
                           biomarker_states, demography, biomarker_map, model_pars, ...){
  return(1)
}

#' Immunity Model For Vaccination Events Only
#' 
#' @description This immunity model should only be used if all exposures are vaccination events. The probability of successful exposure(vaccination event) depends on the number of vaccines an individual has received prior to time t. If the individual is under the maximum vaccinations allotted then the probability of successful exposure event is 1.
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure event
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs 
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure event and biomarker 
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_vacc_events A vector of the maximum number of vaccination events possible for each exposure type; If an exposure type is not a vaccination event then input NA
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure event; If an exposure event is not a vaccination event then input NA
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_only <- function(i, t, x, exposure_histories, 
                              biomarker_states, demography, biomarker_map, 
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
#' @description This immunity model should only be used with vaccines and natural infection. The probability of successful exposure for vaccination events depends on the number of vaccines an individual has received prior to time t and their current age. If the individual is under the maximum vaccinations allotted and is of an age eligible for vaccination then the probability of a successful exposure event is 1. The probability of a successful natural infection event is dependent on the total number of infections that individual has experienced thus far. If the individual is under the maximum number of infections allotted then the probability of a successful exposure event is 1. 
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure events
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs 
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure events and biomarker 
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_events A vector of the maximum number of successful exposure events possible for each exposure event
#' @param vacc_exposures A vector of exposure IDs (x) which represents vaccination events
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure event; If an exposure event is not a vaccination event then input NA
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_ifxn_simple <- function(i, t, x, exposure_histories, 
                                     biomarker_states, demography, biomarker_map, model_pars, max_events, vacc_exposures, vacc_age, ...){
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
    ## Count the total number of successful exposures to x thus far 
    curr_ifx_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
    if(curr_ifx_events<max_events[x]){
      return(1)
    }else{
      return(0)
  }
  }
}



#' Immunity Model For Natural Infection Events With Biomarker Quantity Mediated Protection
#' 
#' @description This immunity model should only be used for natural infection events. The probability of a successful exposure event is dependent on the individual’s biomarker quantity at the time of exposure. User specified biomakrer_prot_midpoint and biomarker_prot_width within model_pars is used to calculate biomarker-mediated protection.   
#'    
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure events
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure events and biomarker
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_ifxn_biomarker_prot <- function(i, t, x, exposure_histories, 
                              biomarker_states, demography, biomarker_map, model_pars, ...){
    ## Find biomarkers which are boosted by this exposure type
    ## The assumption here is that the biomarker quantity will determine if an individual is protected
    b<-biomarker_map$biomarker_id[biomarker_map$exposure_id==x]
    ## Find current biomarker quantity to all relevant biomarkers
    curr_t <- biomarker_states[i,t,b] 
    curr_t <- sum(curr_t)
    
    ## Pull out necessary variables 
    biomarker_prot_midpoint <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="biomarker_prot_midpoint"] 
    biomarker_prot_width <- model_pars$mean[model_pars$exposure_id==x  & model_pars$name=="biomarker_prot_width"]
    
    ## Create a function to calculate the risk of infection at a given biomarker quantity
    biomarker_protection <- function(biomarker_quantity, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(biomarker_quantity - alpha1)))
      return(risk)
    }
    
    prob_success<- (1-biomarker_protection(curr_t, biomarker_prot_midpoint, biomarker_prot_width))
    
    return(prob_success)
  }


#' Immunity Model For Vaccination Events and Natural Infection Events With Biomarker Quantity Mediated Protection
#' 
#' @description This immunity model should be used if exposures represent vaccination and natural infection events. The probability of a successful vaccination exposure event depends on the number of vaccines received prior to time t while the probability of successful infection is dependent on the biomarker quantity at the time of exposure and the total number of successful infections prior to that point.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure events
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param biomarker_map A table specifying the relationship between exposure events and biomarkers
#' @param model_pars A tibble of parameters needed for the immunity model
#' @param max_events A vector of the maximum number of successful exposure events possible for each exposure event; If an exposure type is not a vaccination event then input NA
#' @param vacc_exposures A vector of exposure IDs (x) which represents vaccination events
#' @param vacc_age A vector of the minimum age at which an individual is eligible for vaccination for each exposure event; If an exposure event is not a vaccination event then input NA
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_ifxn_biomarker_prot <- function(i, t, x, exposure_histories,
                           biomarker_states, demography, biomarker_map, model_pars, max_events, vacc_exposures, vacc_age=1, ...){
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
  else {
    ## If the exposure event is an infection event
    ## Count the total number of successful exposures to x thus far 
    curr_ifx_events<-sum(exposure_histories[i,1:t-1,x], na.rm=TRUE)
    
    ## If the current number of successful exposures to x is less than the maximum number of successful exposures
    if(curr_ifx_events<max_events[x]){
      
        ## Find biomarkers which are boosted by this exposure type
        ## The assumption here is that the biomaker quantity will determine if an individual is protected
         b<-biomarker_map$biomarker_id[biomarker_map$exposure_id==x]
        ## Find current biomarker quantity to all relevant biomarkers
        curr_t <- biomarker_states[i,t,b]
        curr_t <- sum(curr_t)
    
        ## Pull out necessary variables 
       biomarker_prot_midpoint <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="biomarker_prot_midpoint"]
        biomarker_prot_width <- model_pars$mean[model_pars$exposure_id==x & model_pars$name=="biomarker_prot_width"]
    
        ## Create a function to calculate the risk of infection at a given biomarker
        biomarker_protection <- function(biomarker_quantity, alpha1, beta1){
          risk <- 1 - 1/(1 + exp(beta1*(biomarker_quantity - alpha1)))
         return(risk)
       }
     
          prob_success<- (1-biomarker_protection(curr_t, biomarker_prot_midpoint, biomarker_prot_width))
    
       return(prob_success)
    }
    else{
      return(0)
    }
  }
}

