#' Here I will go through main_test_JAH.R and fill in the generic calls with
#' more specific models.


## Specify the number of individuals in the simulation 
N <- 10
## Specify the number of time periods to simulate 
times <- seq(1,12,by=1) 
## Specify the number of exposure types 
N_exposure_ids <- 3
## Specify which antigens are present in each exposure type
antigen_map <- tibble(exposure_id=c(1,2,3,3),antigen_id=c(1,2,1,2)) 

## Pre load the demography categories, values and distributions 
aux <- list("SES"=list("name"="SES","options"=c("low","medium","high"), "distribution"=c(0.2,0.2,0.6)),
            "NS"=list("name"="NS","options"=c("low","high"),"distribution"=c(0.5,0.5)),
            "Sex"=list("name"="Sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),
            "location"=list("name"="Location","options"=c("North", "South", "East", "West"), "distribution"=c(0.25,0.25,0.25,0.25)) )

## Simulate demography settings 
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=12, prob_removal=0.9, aux=aux)

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))
## Set observation time
observation_times <- NULL
## Specify the force of infection array
# lambdas <- array(rep(0.01,length(times)), dim=c(length(times),1,1)) ## This doesn't match the FOI structure that was planned 
## Specify the force of infection array (different format)
location<-c("North", "South", "East", "West")
lambdas <- array(rep(0.01,length(location)), dim=c(length(location),max(times),N_exposure_ids))

##Specify the antibody kinetics parameters
theta <- list("boost_mean"=2,"boost_sd"=1)


## Define the exposure model which which will determine if an individual is exposed at a specific location in time.
exposure_model <- function(i, t, e, lambdas, demography){
  l <- demography$Location[demography$i==i & demography$times==t]
  p <- lambdas[l, t, e]
  p
}


## Define the immunity model which will determine if an individual becomes infected if an exposure occured
immunity_model <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map, lambdas, theta,...){
  ## Only look at infection exposures 
  if(antigen_map$class[antigen_map$exposure_id==e]==1){
    ag<-antigen_map$antigen_id[antigen_map$exposure_id==e]
    curr_t <- antibody_states[i,t,ag] ## Find current titre
    
    ## Create a function to calculate the risk of infection at a given titre
    titre_protection <- function(titre, alpha1, beta1){
      risk <- 1 - 1/(1 + exp(beta1*(titre - alpha1)))
      return(risk)
    }
    
    ## Create a function to generate the probability of infection at a given titre
    p_infection <- function(phi, titre, alpha1, beta1){
      p <- phi*(1-titre_protection(titre, alpha1 , beta1))
      return(p)
    }
    
    ## Pull out necessary variables 
    l <- demography$Location[demography$i==i & demography$times==t]
    titre_prot_midpoint <- theta$mean[theta$exposure_id==e & theta$ag==ag & theta$name=="titre_prot_midpoint"]
    titre_prot_width <- theta$mean[theta$exposure_id==e & theta$ag==ag & theta$name=="titre_prot_width"]
    
    ## Find probability of infection given current titre level and protective effect of titres
    prob_infection <- p_infection(lambdas[l,t,e], curr_t, titre_prot_midpoint, titre_prot_width)
    infected <- as.integer(runif(1) < prob_infection)
    return(infected)
    
    
  ## If the exposure event is a vaccination then titre-mediated protection doesn't matter
  } 
  if(antigen_map$class[antigen_map$exposure_id==e]==2){
    return(0)
  }
}

observation_model <- NULL

  

