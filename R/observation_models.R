#' Observation Model Version 1
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and no added noise.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param demography Demography information 
#' @param boundary A list containing the lower bound and upper bound of the assay
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_V1<-function(antibody_states,theta,demography,boundary){
  lower_bound<-boundary[1]
  upper_bound<-boundary[2]
  antibody_states$observed<-antibody_states$value
  antibody_states$observed<-ifelse(antibody_states$observed<lower_bound,0,antibody_states$observed)
  antibody_states$observed<-ifelse(antibody_states$observed>upper_bound,upper_bound,antibody_states$observed)
  antibody_states
}

#' Observation Model Version 2
#' 
#' @description This observation model observes the latent titer values given a discrete assay with user-specified ranges within discrete and no added noise.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param demography Demography information 
#' @param discrete The cutoffs for the discrete assay in a list format
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_V2<-function(antibody_states,theta,demography,discrete){
  discrete_tmp<-c(discrete,Inf)
  antibody_states$observed<-cut(antibody_states$value, breaks=discrete_tmp, right=FALSE, labels=discrete)
  antibody_states
}

#' Observation Model Version 3
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution is defined within theta as the “obs_sd” parameter.
#' 
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param demography Demography information 
#' @param boundary A list containing the lower bound and upper bound of the assay
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_V3<-function(antibody_states,theta,demography,boundary){
  ag_tmp<-unique(antibody_states$ag)
  antibody_states_new<-NULL
  lower_bound<-boundary[1]
  upper_bound<-boundary[2]
  for(ags in seq_along(ag_tmp)){
    if(theta$distribution[theta$antigen_id==ag & theta$name=="obs_sd"]=="log-normal"){
      antibody_states_tmp<-antibody_states %>% filter(ag==ags)
      antibody_states_tmp$observed<-rlnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,theta$mean[theta$antigen_id==ags & theta$name=="obs_sd"])
    }
    if(theta$distribution[theta$antigen_id==ag & theta$name=="obs_sd"]=="normal"){
      antibody_states_tmp<-antibody_states %>% filter(ag==ags)
      antibody_states_tmp$observed<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,theta$mean[theta$antigen_id==ags & theta$name=="obs_sd"])
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  antibody_states<- antibody_states_new %>% arrange(i, t, ag)
  antibody_states$observed<-ifelse(antibody_states$observed<lower_bound,0,antibody_states$observed)
  antibody_states$observed<-ifelse(antibody_states$observed>upper_bound,upper_bound,antibody_states$observed)
  antibody_states$observed<-ifelse(antibody_states$observed<0,0,antibody_states$observed)
  antibody_states
}


#' Observation Model Version 4
#' 
#' @description This observation model observes the latent titer values given a discrete assay with user-specified ranges within discrete and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution is defined within theta as the “obs_sd” parameter.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param demography Demography information 
#' @param discrete The cutoffs for the discrete assay in a list format
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_V4<-function(antibody_states,theta,demography, discrete){
  ag_tmp<-unique(antibody_states$ag)
  discrete_tmp<-c(discrete,Inf)
  antibody_states_new<-NULL
  for(ags in seq_along(ag_tmp)){
    if(theta$distribution[theta$antigen_id==ag & theta$name=="obs_sd"]=="log-normal"){
      antibody_states_tmp<-antibody_states %>% filter(ag==ags)
      antibody_states_tmp$temp<-rlnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,theta$mean[theta$antigen_id==ags & theta$name=="obs_sd"])
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=discrete_tmp, right=FALSE, labels=discrete)
      antibody_states_tmp$temp<-NULL
    }
    if(theta$distribution[theta$antigen_id==ag & theta$name=="obs_sd"]=="normal"){
      antibody_states_tmp<-antibody_states %>% filter(ag==ags)
      antibody_states_tmp$temp<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,theta$mean[theta$antigen_id==ags & theta$name=="obs_sd"])
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=discrete_tmp, right=FALSE, labels=discrete)
      antibody_states_tmp$temp<-NULL
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  antibody_states<- antibody_states_new %>% arrange(i, t, ag)
  antibody_states
}
