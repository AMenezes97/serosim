#' Observation Model For Continuous Assays With No Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with no added noise.Therefore the true titer is the observed titer.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous<-function(antibody_states,model_pars, ...){
  antibody_states$observed<-antibody_states$value
  antibody_states_new
}

#' Observation Model For Continuous Assays With Detection Limits And No Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and no added noise.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param bounds A tibble containing the assay lower bound and upper bound for all antigens; column names=antigen_id, name, and value where name is either "lower_bound" or "upper_bound"
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_bounded<-function(antibody_states,model_pars, bounds, ...){
  antibody_states$observed<-antibody_states$value
  antibody_states_new<-NULL
  for(ags in unique(antibody_states$ag)){ ## For each antigen
    ## Pull out lower and upper bound for the assays
    lower_bound<-bounds$value[bounds$antigen_id==ags & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$antigen_id==ags & bounds$name=="upper_bound"]
    ## Create a new data table containing only the particular antigen
    antibody_states<-data.table(antibody_states)
    antibody_states_tmp<-antibody_states[antibody_states$ag==ags,]
    ## Adjust the observed values given the assay upper and lower bounds
    antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<lower_bound,0,antibody_states_tmp$observed)
    antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed>upper_bound,upper_bound,antibody_states_tmp$observed)
    antibody_states_new<- rbind(antibody_states_new, antibody_states_tmp)
  }
  antibody_states_new
}

#' Observation Model For Discrete Assays With No Added Noise
#' 
#' @description This observation model observes the latent titer values given a discrete assay with user-specified ranges within discrete and no added noise.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each antigen. Each row contains all of the cutoffs for that antigen starting with 0. For example, all of the cutoffs for the assay measuring antigen 1 specific antibodies are in the first row of this matrix.
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_discrete<-function(antibody_states,model_pars, cutoffs, ...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(ags in unique(antibody_states$ag)){ ## For each antigen
    ## Pull out the assay cutoffs 
    cutoffs_ag<-cutoffs[ags,]
    cutoffs_tmp<-c(cutoffs_ag,Inf)
    antibody_states_tmp<-antibody_states[antibody_states$ag==ags,]
    antibody_states_tmp$observed<-cut(antibody_states_tmp$value, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_ag)
    antibody_states_new<- rbind(antibody_states_new, antibody_states_tmp)
  }
  antibody_states_new
}

#' Observation Model For Continuous Assays With Detection Limits And Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter.
#' 
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param bounds A tibble containing the assay lower bound and upper bound for all antigens; column names=antigen_id, name, and value where name is either "lower_bound" or "upper_bound"
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_bounded_noise<-function(antibody_states,model_pars, bounds, ...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(ags in unique(antibody_states$ag)){
    lower_bound<-bounds$value[bounds$antigen_id==ags & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$antigen_id==ags & bounds$name=="upper_bound"]
    antibody_states_tmp<-antibody_states[antibody_states$ag==ags,]
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="log-normal"){
      antibody_states_tmp$observed<-rlnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<lower_bound,0,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed>upper_bound,upper_bound,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$observed<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<lower_bound,0,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed>upper_bound,upper_bound,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, ag)
  observed_states
}

#' Observation Model For Continuous Assays With  Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter.
#' 
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_noise<-function(antibody_states,model_pars, ...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(ags in unique(antibody_states$ag)){
    antibody_states_tmp<-antibody_states[antibody_states$ag==ags,]
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="log-normal"){
      antibody_states_tmp$observed<-rlnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$observed<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, ag)
  observed_states
}

#' Observation Model For Discrete Assays With Added Noise
#' 
#' @description This observation model observes the latent titer values given a discrete assay with user-specified ranges within discrete and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and antigens  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each antigen. Each row contains all of the cutoffs for that antigen starting with 0. For example, all of the cutoffs for the assay measuring antigen 1 specific antibodies are in the first row of this matrix.
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_discrete_noise<-function(antibody_states,model_pars, cutoffs, ...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(ags in unique(antibody_states$ag)){
    cutoffs_ag<-cutoffs[ags,]
    cutoffs_tmp<-c(cutoffs_ag,Inf)
    antibody_states_tmp<-antibody_states[antibody_states$ag==ags,]
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="log-normal"){
      antibody_states_tmp$temp<-rlnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_ag)
      antibody_states_tmp$temp<-NULL
    }
    if(model_pars$distribution[model_pars$antigen_id==ags & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$temp<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$antigen_id==ags & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_ag)
      antibody_states_tmp$temp<-NULL
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, ag)
  observed_states
}

