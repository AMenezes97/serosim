#' Here I will go through main_test_JAH.R and fill in the generic calls with
#' more specific models.
library(serosim)



##******************Component 1: Simulation Settings**************************** 
## Specify the number of time periods to simulate 
times <- seq(1,12,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))



##******************Component 2: Population Demography************************** 

## Specify the number of individuals in the simulation 
N <- 10

## Pre load the demography categories, values and distributions 
aux <- list("SES"=list("name"="SES","options"=c("low","medium","high"), "distribution"=c(0.2,0.2,0.6)),
            "NS"=list("name"="NS","options"=c("low","high"),"distribution"=c(0.5,0.5)),
            "Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),
            "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "distribution"=c(0.25,0.25,0.25,0.25)) )

## Simulate demography settings 
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=12, prob_removal=0.9, aux=aux)




##*****************Component 3: Force of Infection and Exposure Model***********
## Specify the number of exposures (eventually, this won't be needed because the serosim function can extract this from the antigen_map)
N_exposure_ids <- 2

## Specify the force of infection array (different format)
group<-c(1,2,3,4)
lambdas <- array(rep(0.01,length(group)), dim=c(length(group),max(times),N_exposure_ids))


## Define the exposure model which which will determine if an individual is exposed at a specific location in time.
## This is a simple exposure model depending only on the force of infection at that time for that group
exposure_model <- function(i, t, e, lambdas, demography){
  g <- demography$group[demography$i==i & demography$times==t]
  p <- lambdas[g, t, e]
  p
}




##***********************Component 4: Antigen Map*******************************
## Specify which antigens are present in each exposure type
antigen_map <- tibble(exposure_id=c(1,2,3,3),antigen_id=c(1,2,1,2)) 




##**************Optional Component: Antibody Kinetics Parameters****************





##**********************Component 5: Immunity Model*****************************





## Set observation time
observation_times <- NULL


##Specify the antibody kinetics parameters
theta <- list("boost_mean"=2,"boost_sd"=1)


## Specify the max number of vaccination events allowed for each vaccine exposure type
max_vacc_events<- list("1"=3,"2"=1)

## Specify which Exposure_IDs correspond to vaccination events  
vacc_exposures<-c(1,2)

## Define the immunity model which will determine if an exposure is successful 
immunity_model <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map, max_vacc_events, vacc_exposures, lambdas, theta, ...){
  ## If vaccination event (ie. e==vaccination events), then guaranteed exposure unless the individual has already been vaccinated
  if(e %in% c(vacc_exposures)){  	  
    ## Count the total number of successful exposures to e thus far 
    curr_vacc_events<-sum(exposure_histories[i,1:t-1,e])
    ## If number of successful exposures is less than the max number of vaccination events then vaccine exposure is successful 
    if(curr_vacc_events<max_vacc_events[e]){
      return(1)
    }else{
      return(0)
    }
  } else {
    ## Find antigens which are boosted by this exposure type
    ## The assumption here is that the titer levels to these antigens will determine if an individual is protected
    ag<-antigen_map$antigen_id[antigen_map$exposure_id==e]
    ## Find current titer to all relevant antigens
    curr_t <- antibody_states[i,t,ag] ## How to deal with titers against multiple antigens? Should they be added?
    
    ## Pull out necessary variables 
    titer_prot_midpoint <- theta$mean[theta$exposure_id==e & theta$ag==ag & theta$name=="titer_prot_midpoint"]
    titer_prot_width <- theta$mean[theta$exposure_id==e & theta$ag==ag & theta$name=="titer_prot_width"]
    
    ## Create a function to calculate the risk of infection at a given titer
    titer_protection <- function(titer, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(titer - alpha1)))
      return(risk)
    }
    
    prob_success<- (1-titer_protection(curr_t, titer_prot_midpoint, titer_prot_width))
    
    return(prob_success)
  }
}

  
    
    
    
observation_model <- function(antibody_states, theta, demography, ...){
  ## antibody_states needs to be converted to the long format first
  collapse_array <- function(titers, N_antigen_ids){
    collapsed<-NULL
    for(antigen in 1:N_antigen_ids){
      tmp <- reshape2::melt(titers[,,antigen])
      colnames(tmp) <- c("Individual","Time","Titer")
      tmp$Antigen <- antigen
      collapsed<-rbind(collapsed,tmp)
    }
    return(collapsed)
  }
  
  antibody_states$observed <- rnorm(nrow(antibody_states),antibody_states$value,theta[["obs_sd"]])
  antibody_states
}



draw_parameters <- function(i, t, e, ag, demography, theta, antibody_state, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    if(theta_tmp$distribution == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
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
      
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(theta_tmp$mean[par], theta_tmp$sd[par]), normal_to_lognormal_sd(theta_tmp$mean[par], theta_tmp$sd[par]))
      par_names[par] <- theta_tmp$name[par]
    }
    else{
      pars[par] <- rnorm(1, theta_tmp$mean[par], theta_tmp$sd[par])
      par_names[par] <- theta_tmp$name[par]
    }
    #Add titer-dependent boosting  
    if(par_names[par] %in% c("boost_short","boost_long")){
      titer_threshold <- min(antibody_states[i,t1,ag], theta_tmp[theta_tmp$name=="titer_ceiling_threshold" & theta_tmp$ag==ag, "mean"])
      pars[par] <- pars[par]*(1-theta_tmp[theta_tmp$name=="titer_ceiling_gradient" & theta_tmp$ag==ag, "mean"]*titer_threshold)
    }
  }
  all_pars <- tibble(i=i, t=t, e=e, ag=ag, name=par_names, value=pars)
  return(all_pars)
}
  

