
devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)

## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))

N<-500

demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=120, prob_removal=0)

antigen_map <- tibble(exposure_id=c(1,2),antigen_id=c(1,1)) 

## Create an empty array to store the force of infection for all exposure types
lambdas <- array(0, dim=c(1,max(times),n_distinct(antigen_map$exposure_id)))

## Specify the force of infection for exposure ID 1 which represents natural infection
lambdas[,,1] <- 0.2

## Specify the force of vaccination for exposure ID 2 which represents vaccination
## Assume this is constant across all time steps
lambdas[,,2] <- 0.4

## Specify a simple exposure model which calculates the probability of exposure directly from the force of infection at that time step
## The probability of exposure is 1-exp(-FOI) where FOI is the force of infection at that time
exposure_model<-exposure_model_simple_FOI

## Specify immunity model within serosim function below 
immunity_model<-immunity_model_vacc_ifxn_titer_prot

## Specify which exposure IDs represent vaccination events 
vacc_exposures<-2

## Specify the age at which an individual is eligible for vaccination (9 months old for measles)
vacc_age<-9

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
max_vacc_events<-c(NA,1)

## Specify the antibody model 
antibody_model<-antibody_model_biphasic

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model, the titer-ceiling parameters needed for draw_parameters and the observation error needed for the observation model are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
theta_path <- system.file("extdata", "theta_cs1.csv", package = "serosim")
theta <- read.csv(file = theta_path, header = TRUE)

## Specify the draw_parameters function to use 
draw_parameters<-draw_parameters_random_fx_titer_dep


## Limits of detection for continuous assays
boundary<-c(2,20)

## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify observation_times to observe antigen 1 (aka Measles antibody titer) across all individuals at the midpoint and the end of the simulation (t=60 and t=120)
obs1 <- tibble(i=1:N,t=60, ag=1)
obs2 <- tibble(i=1:N,t=120, ag=1)
observation_times<-rbind(obs1,obs2)


res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  lambdas,
  antigen_map,
  theta,
  exposure_model,
  immunity_model,
  antibody_model,
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_random_fx_titer_dep,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
)

plot_titers(res$antibody_states)
plot_exposure_prob(res$exposure_probabilities_long)
plot_obs_titers_one_sample(res$observed_antibody_states)
plot_obs_titers_paired_sample(res$observed_antibody_states)
plot_exposure_histories(res$exposure_histories_long)
