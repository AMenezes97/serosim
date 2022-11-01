## The following script will be used to generate example arguments 
## to be stored within serosim and used for examples within the help function.
## This script will follow the same format as README.Rmd
## After the simulation, runserosim outputs will also be stored.
devtools::document("~/Documents/GitHub/serosim")
devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)



## Specify the number of time steps in the simulation
times <- seq(1,120,by=1) 
## Set simulation settings argument needed for runserosim 
simulation_settings <- list("t_start"=1,"t_end"=max(times))



## Create an example demography tibble
example_demography <- generate_pop_demography(N=100, times=times, removal_min=0, removal_max=120, prob_removal=0)

## Create an example antigen map for 2 exposures with the same antigen
example_antigen_map <- tibble(exposure_id=c(1,2),antigen_id=c(1,1)) 

## Create an example foe_pars array
foe_pars <- array(0, dim=c(1,max(times),n_distinct(example_antigen_map$exposure_id)))
foe_pars[,,1] <- 0.2
foe_pars[,,2] <- 0.4
example_foe_pars <-foe_pars



## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model (Section 1.5), the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model (Section 1.7) are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
theta_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
example_model_pars <- read.csv(file = theta_path, header = TRUE)



## Specify observation_times (serological survey sampling design) to observe antigen 1 across all individuals at the end of the simulation (t=120)
observation_times<- tibble(i=1:max(example_demography$i),t=120, ag=1)


res<- runserosim(
  simulation_settings=simulation_settings,
  demography=example_demography,
  observation_times=observation_times,
  foe_pars=example_foe_pars,
  antigen_map=example_antigen_map,
  theta=example_model_pars,
  exposure_model=exposure_model_simple_FOE,
  immunity_model=immunity_model_vacc_ifxn_simple,
  antibody_model=antibody_model_monophasic,
  observation_model=observation_model_continuous_noise,
  draw_parameters=draw_parameters_random_fx,
  
  ## Other arguments needed
  max_events=c(1,1),
  vacc_exposures=2,
  vacc_age=c(NA,9),
)


## Save select runserosim outputs 
example_exposure_histories<-res$exposure_histories_long
example_exposure_probabilities<-res$exposure_probabilities_long
example_antibody_states<-res$antibody_states
example_observed_antibody_states<-res$observed_antibody_states




## Save all generated data
save(example_demography, 
     example_antigen_map,
     example_foe_pars,
     example_model_pars,
     example_exposure_histories,
     example_exposure_probabilities,
     example_antibody_states,
     example_observed_antibody_states,
     file="~/Documents/GitHub/serosim/data/example.RData")


