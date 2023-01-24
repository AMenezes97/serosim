#' Observation model for continuous assays with no added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with no added noise. Therefore the observed biomarker quantity is simply given by the true latent biomarker quantity.
#'
#' @param biomarker_states tibble containing true biomarker quantities for all individuals across all time steps and biomarkers. Variables should include: 1) i: the individual ID; 2) t: the time period; 3) b: the biomarker ID; 4) value: the latent biomarker quantity for the given i, t and b
#' @param model_pars a tibble containing information for all parameters needed to simulate the observation process. This should usually contain: 1) exposure_id: numeric exposure ID; 2) biomarker_id: numeric biomarker ID; 3) name: the character name of the parameter; 4) mean: numeric mean of this parameter distribution; 5) sd: the numeric standard deviation of the parameter distribution
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @family observation_model
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
#' @family observation_model
#' @export
#'
#' @examples
#' bounds <- tibble(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(2,8))
#' observation_model_continuous_bounded(example_biomarker_states, NULL,bounds)
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

#' Observation model for discrete assays with no added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a discrete assay with user-specified ranges within the "discrete" argument and no added noise.
#'
#' @inheritParams observation_model_continuous
#' @param cutoffs a matrix containing the assay cutoffs for each biomarker. Each row contains all of the cutoffs for that biomarker starting with 0. For example, all of the cutoffs for the assay measuring biomarker 1 quantity are in the first row of this matrix.
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @family observation_model
#' @export
#'
#' @examples
#' breaks <- seq(0,8,by=1)
#' cutoffs <- matrix(breaks,nrow=1,ncol=length(breaks))
#' observation_model_discrete(example_biomarker_states, NULL,cutoffs)
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

#' Observation model for continuous assays with detection limits and added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with user-specified lower and upper limits and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within `model_pars` as the `obs_sd` parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. False negatives are simulated by setting an observed quantity to the assay's lower bound with probability `sensitivity`. False positives are simulated by drawing a random quantity from the bounded range for a true 0 biomarker quantity with probability 1-`specificity`.
#' 
#' @inheritParams observation_model_continuous_bounded
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @family observation_model
#' @export
#'
#' @examples
#' bounds <- tibble(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(2,8))
#' observation_model_continuous_bounded_noise(example_biomarker_states, example_model_pars_numeric, bounds,0.95,0.99)
observation_model_continuous_bounded_noise<-function(biomarker_states,model_pars, bounds, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
    lower_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="lower_bound"]
    upper_bound<-bounds$value[bounds$biomarker_id==bs & bounds$name=="upper_bound"]
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      biomarker_states_tmp$observed<-rlnorm(nrow(biomarker_states_tmp),normal_to_lognormal_mean(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<lower_bound,lower_bound,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed>upper_bound,upper_bound,biomarker_states_tmp$observed)
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      biomarker_states_tmp$observed<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed<lower_bound,lower_bound,biomarker_states_tmp$observed)
      biomarker_states_tmp$observed<-ifelse(biomarker_states_tmp$observed>upper_bound,upper_bound,biomarker_states_tmp$observed)
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed biomarker quantities 
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(biomarker_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
            biomarker_states_tmp$observed[indiv]<-runif(1,lower_bound,upper_bound)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(biomarker_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            biomarker_states_tmp$observed[indiv]<-lower_bound #individual is now incorrectly labeled as a false negative
          }
        }
      }
    }
    biomarker_states_new<-rbind(biomarker_states_new,biomarker_states_tmp)
  }
  observed_states<- biomarker_states_new %>% arrange(i, t, b)
  observed_states
}

#' Observation model for continuous assays with no detection limits and added noise
#' 
#' @description This observation model observes the latent biomarker quantities given a continuous assay with added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the `obs_sd` parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. False negatives are simulated by setting an observed quantity to 0 with probability `sensitivity`. False positives are simulated by drawing a random quantity from the bounded range for a true 0 biomarker quantity with probability 1-`specificity`.
#' 
#' @inheritParams observation_model_continuous_noise
#' @param sensitivity number between 0 and 1 to describe the assay's sensitivity; defaults to 1
#' @param specificity number between 0 and 1 to describe the assay's specificity; defaults to 1
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @family observation_model
#' @export
#'
#' @examples
#' observation_model_continuous_noise(example_biomarker_states, example_model_pars_numeric, 0.95,0.99)
observation_model_continuous_noise<-function(biomarker_states,model_pars, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
    biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      biomarker_states_tmp$observed<-rlnorm(nrow(biomarker_states_tmp),normal_to_lognormal_mean(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
    }
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
      biomarker_states_tmp$observed<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
    }
    if(sensitivity!=1 | specificity!=1){
      ## Calculate the range of observed biomarker quantities
        lower_bound <- min(biomarker_states_tmp$value,na.rm=TRUE)
        upper_bound <- max(biomarker_states_tmp$value,na.rm=TRUE)
        
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
      ## If the individual is truly seronegative take into account the assay specificity 
      if(biomarker_states_tmp$value[indiv]==0){
        if(runif(1)>specificity){
            biomarker_states_tmp$observed[indiv]<-runif(1,lower_bound,upper_bound)
        }
      }
      ## If the individual is truly seropositve take into account the assay sensitivity 
      if(biomarker_states_tmp$value[indiv]!=0){
        if(runif(1)>sensitivity){
          biomarker_states_tmp$observed[indiv]<-0 #individual is now incorrectly labeled as a false negative
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
#' @description  This observation model observes the latent biomarker quantities given a discrete assay with user-specified ranges within discrete and added noise. The added noise represents assay variability and is done by sampling from a distribution with the latent biomarker quantity as the mean and the measurement error as the standard deviation. The observation standard deviation and distribution are defined within model_pars as the “obs_sd” parameter. The user can also use the optional sensitivity and specificity arguments to account for assay sensitivity and specificity. False negatives are simulated by setting an observed quantity to the assay's lower bound with probability `sensitivity`. False positives are simulated by drawing a random quantity from the bounded range for a true 0 biomarker quantity with probability 1-`specificity`.
#'
#' @inheritParams observation_model_continuous_noise
#' @inheritParams observation_model_discrete
#' @param ... 
#'
#' @return `biomarker_states` is returned with a new column, `observed`, for observed biomarker quantities
#' @family observation_model
#' @export
#'
#' @examples
#' breaks <- seq(0,8,by=1)
#' cutoffs <- matrix(breaks,nrow=1,ncol=length(breaks))
#' tmp_pars <- example_model_pars_numeric %>% mutate(sd=ifelse(name=="obs_sd",2,sd))
#' observation_model_discrete_noise(example_biomarker_states, tmp_pars, cutoffs, 0.95,0.99)
observation_model_discrete_noise<-function(biomarker_states,model_pars, cutoffs, sensitivity=1, specificity=1,...){
  biomarker_states_new<-NULL
  biomarker_states<-data.table(biomarker_states)
  for(bs in unique(biomarker_states$b)){
      cutoffs_b<-cutoffs[bs,]
      cutoffs_tmp<-c(cutoffs_b,Inf)
      biomarker_states_tmp<-biomarker_states[biomarker_states$b==bs,]
      
      if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="normal"){
          biomarker_states_tmp$temp<-rnorm(nrow(biomarker_states_tmp),biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"])
      }
      
    if(model_pars$distribution[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]=="log-normal"){
      biomarker_states_tmp$temp<-rlnorm(nrow(biomarker_states_tmp),normal_to_lognormal_mean(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]),normal_to_lognormal_sd(biomarker_states_tmp$value,model_pars$sd[model_pars$biomarker_id==bs & model_pars$name=="obs_sd"]))
    }
      min_observable <- min(as.numeric(cutoffs_b),na.rm=TRUE)
      biomarker_states_tmp$temp <- ifelse( biomarker_states_tmp$temp  < min_observable,min_observable,biomarker_states_tmp$temp)
      biomarker_states_tmp$observed<-cut(biomarker_states_tmp$temp, breaks=cutoffs_tmp, right=FALSE, labels=cutoffs_b)
      biomarker_states_tmp$temp<-NULL
    
    if(sensitivity!=1 | specificity!=1){
      ## Step through each individual 
      for(indiv in 1:length(biomarker_states_tmp$value)){
        ## If the individual is truly seronegative take into account the assay specificity 
        if(biomarker_states_tmp$value[indiv]==0){
          if(runif(1)>specificity){
              biomarker_states_tmp$observed[indiv]<-sample(cutoffs_b,1)
          }
        }
        ## If the individual is truly seropositve take into account the assay sensitivity 
        if(biomarker_states_tmp$value[indiv]!=0){
          if(runif(1)>sensitivity){
            biomarker_states_tmp$observed[indiv]<-cutoffs_b[1] #individual is now incorrectly labeled as having the lowest possible biomarker quantity (false negative)
          }
        }
      }
    }
    biomarker_states_new<-rbind(biomarker_states_new,biomarker_states_tmp)
  }
  observed_states<- biomarker_states_new %>% arrange(i, t, b)
  observed_states
}

