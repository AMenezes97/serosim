#' Observation Model For Continuous Assays With No Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with no added noise.Therefore the true titer is the observed titer.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous<-function(antibody_states,model_pars, ...){
  antibody_states$observed<-antibody_states$value
  return(antibody_states)
}

#' Observation Model For Continuous Assays With Detection Limits And No Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and no added noise.
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param bounds A tibble containing the assay lower bound and upper bound for all biomarkers; column names=biomarker_id, name, and value where name is either "lower_bound" or "upper_bound"
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_bounded<-function(antibody_states,model_pars, bounds, ...){
  antibody_states$observed<-antibody_states$value
  antibody_states_new<-NULL
  for(bs in unique(antibody_states$b)){ ## For each biomarker
    ## Pull out lower and upper bound for the assays
    lower_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="upper_bound"]
    ## Create a new data table containing only the particular biomarker
    antibody_states<-data.table(antibody_states)
    antibody_states_tmp<-antibody_states[antibody_states$b==bs,]
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
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each biomarker. Each row contains all of the cutoffs for that biomarker starting with 0. For example, all of the cutoffs for the assay measuring biomarker 1 specific antibodies are in the first row of this matrix.
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_discrete<-function(antibody_states,model_pars, cutoffs, ...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(bs in unique(antibody_states$b)){ ## For each biomarker
    ## Pull out the assay cutoffs 
    cutoffs_b<-cutoffs[bs,]
    cutoffs_tmp<-c(cutoffs_b,Inf)
    antibody_states_tmp<-antibody_states[antibody_states$b==bs,]
    antibody_states_tmp$observed<-cut(antibody_states_tmp$value, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
    antibody_states_new<- rbind(antibody_states_new, antibody_states_tmp)
  }
  antibody_states_new
}

#' Observation Model For Continuous Assays With Detection Limits And Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with user-specified lower and upper limits and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#' 
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param bounds A tibble containing the assay lower bound and upper bound for all biomarkers; column names=biomarker_id, name, and value where name is either "lower_bound" or "upper_bound"
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_bounded_noise<-function(antibody_states,model_pars, bounds, sensitivity=1, specificity=1,...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(bs in unique(antibody_states$b)){
    lower_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="upper_bound"]
    antibody_states_tmp<-antibody_states[antibody_states$b==bs,]
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      ## Create functions to convert normal distributions to log-normal distributions
      normal_to_lognormal_mean <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        meanlog <- log(normmean ^ 2 / phi)
        return(meanlog)
      }
      
      normal_to_lognormal_sd <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
        return(sdlog)
      }
      
      antibody_states_tmp$observed<-rlnorm(nrow(antibody_states_tmp),normal_to_lognormal_mean(antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<lower_bound,0,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed>upper_bound,upper_bound,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$observed<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<lower_bound,0,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed>upper_bound,upper_bound,antibody_states_tmp$observed)
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed titers for the biomarker
      pos<-antibody_states_tmp %>% filter(antibody_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(antibody_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(antibody_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
            antibody_states_tmp$observed[indiv]<-sample(min:max,size=1)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(antibody_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            antibody_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
          }
        }
      }
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, b)
  observed_states
}

#' Observation Model For Continuous Assays With  Added Noise
#' 
#' @description This observation model observes the latent titer values given a continuous assay with added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter.The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#' 
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_continuous_noise<-function(antibody_states,model_pars, sensitivity=1, specificity=1,...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(bs in unique(antibody_states$b)){
    antibody_states_tmp<-antibody_states[antibody_states$b==bs,]
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      ## Create functions to convert normal distributions to log-normal distributions
      normal_to_lognormal_mean <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        meanlog <- log(normmean ^ 2 / phi)
        return(meanlog)
      }
      
      normal_to_lognormal_sd <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
        return(sdlog)
      }
      
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$observed<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-ifelse(antibody_states_tmp$observed<0,0,antibody_states_tmp$observed)
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed titers for the biomarker
      pos<-antibody_states_tmp %>% filter(antibody_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(antibody_states_tmp$value)){
      ## If the individual is truly seronegative take into account the assay specificity 
      if(antibody_states_tmp$value[indiv]==0){
        if(runif(1)>specificity){
          antibody_states_tmp$observed[indiv]<-sample(min:max,size=1)
        }
      }
      ## If the individual is truly seropositve take into account the assay sensitivity 
      if(antibody_states_tmp$value[indiv]!=0){
        if(runif(1)>sensitivity){
          antibody_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
        }
      }
      }
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, b)
  observed_states
    }

#' Observation Model For Discrete Assays With Added Noise
#' 
#' @description This observation model observes the latent titer values given a discrete assay with user-specified ranges within discrete and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent antibody titer as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter.The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#'
#' @param antibody_states True antibody titers for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each biomarker. Each row contains all of the cutoffs for that biomarker starting with 0. For example, all of the cutoffs for the assay measuring biomarker 1 specific antibodies are in the first row of this matrix.
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return antibody_states is returned with a new column for observed titers
#' @export
#'
#' @examples
observation_model_discrete_noise<-function(antibody_states,model_pars, cutoffs, sensitivity=1, specificity=1,...){
  antibody_states_new<-NULL
  antibody_states<-data.table(antibody_states)
  for(bs in unique(antibody_states$b)){
    cutoffs_b<-cutoffs[bs,]
    cutoffs_tmp<-c(cutoffs_b,Inf)
    antibody_states_tmp<-antibody_states[antibody_states$b==bs,]
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      ## Create functions to convert normal distributions to log-normal distributions
      normal_to_lognormal_mean <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        meanlog <- log(normmean ^ 2 / phi)
        return(meanlog)
      }
      
      normal_to_lognormal_sd <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
        return(sdlog)
      }
      antibody_states_tmp$temp<-rlnorm(nrow(antibody_states_tmp),normal_to_lognormal_mean(antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
      antibody_states_tmp$temp<-NULL
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      antibody_states_tmp$temp<-rnorm(nrow(antibody_states_tmp),antibody_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      antibody_states_tmp$observed<-cut(antibody_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
      antibody_states_tmp$temp<-NULL
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed titers for the biomarker
      pos<-antibody_states_tmp %>% filter(antibody_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(antibody_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(antibody_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
            antibody_states_tmp$observed[indiv]<-sample(min:max,size=1)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(antibody_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            antibody_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
          }
        }
      }
    }
    antibody_states_new<-rbind(antibody_states_new,antibody_states_tmp)
  }
  observed_states<- antibody_states_new %>% arrange(i, t, b)
  observed_states
}

