#' Immunity Model Where All Exposures Are Successful
#'
#' @description Probability of success is 1 so all exposure events are successful
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param demography Demography information 
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param theta Tibble including titer-mediated protection parameters 
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_all_successful <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map, theta, ...){
  return(1)
}

#' Immunity Model For Vaccination Events Only
#' 
#' @description This immunity model should only be used if all exposures are vaccination events. The probability of successful exposure(vaccination event) depends on the number of vaccines an individual has received prior to time t. If the individual is under the maximum vaccinations allotted then the probability of successful exposure is 1.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param demography Demography information 
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param theta Tibble including titer-mediated protection parameters 
#' @param max_vacc_events A list of the maximum number of vaccination events possible for each exposure type
#' @param vacc_age The minimum age at which an individual is eligible for vaccination; defaults to 1 
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_only <- function(i, t, e, exposure_histories, 
                              antibody_states, demography, antigen_map, 
                              theta, max_vacc_events, vacc_age=1,...){
  ## Calculate the individual's current age
  birth_time<-unique(demography$birth[demography$i==i])
  curr_age<- t-birth_time
  
  ## If the individual is above the minimum age of vaccination 
  if(curr_age>=vacc_age){
  ## Count the total number of successful exposures to e thus far 
    curr_vacc_events<-sum(exposure_histories[i,1:t-1,e], na.rm=TRUE)
    ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
    if(curr_vacc_events<max_vacc_events[e]){
      return(1)
    }else{
      return(0)
    }
  }else{
    return(0)
  }
}
 

#' Immunity Model For Vaccination Events and Successful Natural Infection Events
#' 
#' @description This immunity model should only be used with vaccine and natural infection. The probability of successful exposure(vaccination event) depends on the number of vaccines an individual has received prior to time t. If the individual is under the maximum vaccinations allotted then the probability of successful exposure is 1. The probability of a successful natural infection is 1.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param demography Demography information 
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param theta Tibble including titer-mediated protection parameters 
#' @param max_vacc_events A list of the maximum number of vaccination events possible for each exposure type
#' @param vacc_exposures A list of exposure IDs (e) which represents vaccination events
#' @param vacc_age The minimum age at which an individual is eligible for vaccination; defaults to 1 
#' @param ... 
#'
#' @return  A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_successful_ifxn <- function(i, t, e, exposure_histories, 
                                     antibody_states, demography, antigen_map, theta, max_vacc_events, vacc_exposures, vacc_age=1, ...){
  ## If an exposure event is a vaccination event, then guaranteed exposure unless the individual has already been vaccinated
  if(e %in% c(vacc_exposures)){  	
    ## Calculate the individual's current age
    birth_time<-unique(demography$birth[demography$i==i])
    curr_age<- t-birth_time
    
    ## If the individual is above the minimum age of vaccination 
    if(curr_age>=vacc_age){
         ## Count the total number of successful exposures to e thus far 
          curr_vacc_events<-sum(exposure_histories[i,1:t-1,e], na.rm=TRUE)
        ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
        if(curr_vacc_events<max_vacc_events[e]){
           return(1)
             }else{
             return(0)
              }
    }
    if(curr_age<vacc_age){
      return(0)
    }
  }
  else{
    return(1)
  }
}



#' Immunity Model For Natural Infection Events With Titer-Mediated Protection
#' 
#' @description  This immunity model should only be used if all exposures are natural infection events. The probability of successful exposure is dependent on the individual’s antibody titer at the time of exposure. 
#'    
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param demography Demography information 
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param theta Tibble including titer-mediated protection parameters 
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_ifxn_titer_prot <- function(i, t, e, exposure_histories, 
                              antibody_states, demography, antigen_map, theta, ...){
    ## Find antigens which are boosted by this exposure type
    ## The assumption here is that the titer levels to these antigens will determine if an individual is protected
    ag<-antigen_map$antigen_id[antigen_map$exposure_id==e]
    ## Find current titer to all relevant antigens
    curr_t <- antibody_states[i,t,ag] ## How to deal with titers against multiple antigens? Should they be added?
    curr_t <- sum(curr_t)
    
    ## Pull out necessary variables 
    titer_prot_midpoint <- theta$mean[theta$exposure_id==e & theta$name=="titer_prot_midpoint"] ## Would each antigen have it's own value?
    titer_prot_width <- theta$mean[theta$exposure_id==e  & theta$name=="titer_prot_width"]
    
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
#' @param e exposure
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param demography Demography information 
#' @param antigen_map Object determining relationship between exposure IDs and antigens
#' @param theta Tibble including titer-mediated protection parameters 
#' @param max_vacc_events A list of the maximum number of vaccination events possible for each exposure type
#' @param vacc_exposures A list of exposure IDs (e) which represents vaccination events
#' @param vacc_age The minimum age at which an individual is eligible for vaccination; defaults to 1 
#' @param ... 
#'
#' @return A probability of successful exposure is returned
#' @export
#'
#' @examples
immunity_model_vacc_ifxn_titer_prot <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map, theta, max_vacc_events, vacc_exposures, vacc_age=1, ...){
  ## If an exposure event is a vaccination event, then guaranteed exposure unless the individual has already been vaccinated
  if(e %in% c(vacc_exposures)){  
    ## Calculate the individual's current age
    birth_time<-unique(demography$birth[demography$i==i])
    curr_age<- t-birth_time
    
    ## If the individual is above the minimum age of vaccination 
    if(curr_age>=vacc_age){
    ## Count the total number of successful exposures to e thus far 
    curr_vacc_events<-sum(exposure_histories[i,1:t-1,e], na.rm=TRUE)
    ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
    if(curr_vacc_events<max_vacc_events[e]){
      return(1)
    }else{
      return(0)
    }
    if(curr_age<vacc_age){
      return(0)
    }
  } 
    }
  else {
    ## Find antigens which are boosted by this exposure type
    ## The assumption here is that the titer levels to these antigens will determine if an individual is protected
    ag<-antigen_map$antigen_id[antigen_map$exposure_id==e]
    ## Find current titer to all relevant antigens
    curr_t <- antibody_states[i,t,ag] ## How to deal with titers against multiple antigens? Should they be added?
    curr_t <- sum(curr_t)
    
    ## Pull out necessary variables 
    titer_prot_midpoint <- theta$mean[theta$exposure_id==e & theta$name=="titer_prot_midpoint"]
    titer_prot_width <- theta$mean[theta$exposure_id==e & theta$name=="titer_prot_width"]
    
    ## Create a function to calculate the risk of infection at a given titer
    titer_protection <- function(titer, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(titer - alpha1)))
      return(risk)
    }
    
    prob_success<- (1-titer_protection(curr_t, titer_prot_midpoint, titer_prot_width))
    
    return(prob_success)
  }
}

