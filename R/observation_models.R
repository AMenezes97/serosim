#' Observation model for continuous assays with no added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with no added noise. Therefore the observed biomarker quantity is simply given by the true latent biomarker quantity.
#'
#' @param biomarker_states tibble containing true biomarker quantities for all individuals across all time steps and biomarkers. Variables should include: 1) i: the individual ID; 2) t: the time period; 3) b: the biomarker ID; 4) value: the latent biomarker quantity for the given i, t and b
#' @param model_pars a tibble containing information for all parameters needed to simulate the observation process. This should usually contain: 1) exposure_id: numeric exposure ID; 2) biomarker_id: numeric biomarker ID; 3) name: the character name of the parameter; 4) mean: numeric mean of this parameter distribution; 5) sd: the numeric standard deviation of the parameter distribution
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @export
#'
#' @examples
#' observation_model_continuous(example_biomarker_states, NULL)
observation_model_continuous<-function(biomarker_states,model_pars, ...){
  biomarker_states$observed<-biomarker_states$value
  return(biomarker_states)
}

#' Observation model for continuous assays with detection limits and no added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with user-specified lower and upper limits and no added noise.
#'
#' @inheritParams observation_model_continuous
#' @param bounds a tibble containing the assay lower bound and upper bound for all biomarkers; column namesare 1) biomarker_id; 2) name; 3) value, where name is either `lower_bound` or `upper_bound`
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @export
#'
#' @examples
#' bounds <- tibble(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(2,8))
#' observation_model_continuous(example_biomarker_states, NULL,bounds)
observation_model_continuous_bounded<-function(biomarker_states,model_pars, bounds, ...){
  biomarker_states$observed<-biomarker_states$value
  biomarker_states_new<-NULL
  for(bs in unique(biomarker_states$b)){ ## For each biomarker
    ## Pull out lower and upper bound for the assays
    lower_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="upper_bound"]
    ## Create a new data table containing only the particular biomarker
    biomarker_states<-data.table(biomarker_states)
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
    ## Adjust the observed values given the assay upper and lower bounds
    biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<lower_bound,lower_bound,biomarker_states_tmp$observed)
    biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed>upper_bound,upper_bound,biomarker_states_tmp$observed)
    biomarker_states_new<- rbind(biomarker_states_new, biomarker_states_tmp)
  }
  biomarker_states_new
}

#' Observation Model For Discrete Assays With No Added Noise
#' 
#' @description This observation model observes the latent biomarker quantities given a discrete assay with user-specified ranges within "discrete" argument and no added noise.
#'
#' @param biomarker_states True biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each biomarker. Each row contains all of the cutoffs for that biomarker starting with 0. For example, all of the cutoffs for the assay measuring biomarker 1 quantity are in the first row of this matrix.
#' @param ... 
#'
#' @return biomarker_states is returned with a new column for observed biomarker quantities
#' @export
#'
#' @examples
observation_model_discrete<-function(biomarker_states,model_pars, cutoffs, ...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){ ## For each biomarker
    ## Pull out the assay cutoffs 
    cutoffs_b<-cutoffs[bs,]
    cutoffs_tmp<-c(cutoffs_b,Inf)
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
    biomarker_states_tmp$observed<-cut(biomarker_states_tmp$value, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
    biomarker_states_new<- rbind(biomarker_states_new, biomarker_states_tmp)
  }
  biomarker_states_new
}

#' Observation Model For Continuous Assays With Detection Limits And Added Noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with user-specified lower and upper limits and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#' 
#' @param biomarker_states True biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param bounds A tibble containing the assay lower bound and upper bound for all biomarkers; column names=biomarker_id, name, and value where name is either "lower_bound" or "upper_bound"
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return biomarker_states is returned with a new column for observed biomarker quantities
#' @export
#'
#' @examples
observation_model_continuous_bounded_noise<-function(biomarker_states,model_pars, bounds, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
    lower_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="upper_bound"]
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
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
      
      biomarker_states_tmp$observed<-rlnorm(nrow(biomarker_states_tmp),normal_to_lognormal_mean(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<lower_bound,0,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed>upper_bound,upper_bound,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<0,0,biomarker_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      biomarker_states_tmp$observed<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<lower_bound,0,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed>upper_bound,upper_bound,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<0,0,biomarker_states_tmp$observed)
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed biomarker quantities 
      pos<-biomarker_states_tmp %>% filter(biomarker_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(biomarker_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
            biomarker_states_tmp$observed[indiv]<-sample(min:max,size=1)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(biomarker_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            biomarker_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
          }
        }
      }
    }
    biomarker_states_new<-rbind(biomarker_states_new,biomarker_states_tmp)
  }
  observed_states<- biomarker_states_new %>% arrange(i, t, b)
  observed_states
}

#' Observation Model For Continuous Assays With Added Noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#' 
#' @param biomarker_states True biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return biomarker_states is returned with a new column for observed biomarker quantities
#' @export
#'
#' @examples
observation_model_continuous_noise<-function(biomarker_states,model_pars, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
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
      
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<0,0,biomarker_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      biomarker_states_tmp$observed<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<0,0,biomarker_states_tmp$observed)
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed biomarker quantities
      pos<-biomarker_states_tmp %>% filter(biomarker_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
      ## If the individual is truly seronegative take into account the assay specificity 
      if(biomarker_states_tmp$value[indiv]==0){
        if(runif(1)>specificity){
          biomarker_states_tmp$observed[indiv]<-sample(min:max,size=1)
        }
      }
      ## If the individual is truly seropositve take into account the assay sensitivity 
      if(biomarker_states_tmp$value[indiv]!=0){
        if(runif(1)>sensitivity){
          biomarker_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
        }
      }
      }
    }
    biomarker_states_new<-rbind(biomarker_states_new,biomarker_states_tmp)
  }
  observed_states<- biomarker_states_new %>% arrange(i, t, b)
  observed_states
    }

#' Observation Model For Discrete Assays With Added Noise
#' 
#' @description  This observation model observes the latent biomarker quantities given a discrete assay with user-specified ranges within discrete and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. 
#'
#' @param biomarker_states True biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of observation model parameters 
#' @param cutoffs A matrix containing the assay cutoffs for each biomarker. Each row contains all of the cutoffs for that biomarker starting with 0. For example, all of the cutoffs for the assay measuring biomarker 1 quantity are in the first row of this matrix.
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return biomarker_states is returned with a new column for observed biomarker quantities
#' @export
#'
#' @examples
observation_model_discrete_noise<-function(biomarker_states,model_pars, cutoffs, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
    cutoffs_b<-cutoffs[bs,]
    cutoffs_tmp<-c(cutoffs_b,Inf)
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
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
      biomarker_states_tmp$temp<-rlnorm(nrow(biomarker_states_tmp),normal_to_lognormal_mean(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
      biomarker_states_tmp$observed<-cut(biomarker_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
      biomarker_states_tmp$temp<-NULL
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      biomarker_states_tmp$temp<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      biomarker_states_tmp$observed<-cut(biomarker_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
      biomarker_states_tmp$temp<-NULL
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed quantites for the biomarker
      pos<-biomarker_states_tmp %>% filter(biomarker_states_tmp$observed>0)
      pos<-pos$observed
      min<-min(pos)
      max<-max(pos)
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(biomarker_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
            biomarker_states_tmp$observed[indiv]<-sample(min:max,size=1)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(biomarker_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            biomarker_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as seronegative
          }
        }
      }
    }
    biomarker_states_new<-rbind(biomarker_states_new,biomarker_states_tmp)
  }
  observed_states<- biomarker_states_new %>% arrange(i, t, b)
  observed_states
}

