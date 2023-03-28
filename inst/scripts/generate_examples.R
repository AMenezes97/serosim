## The following script will be used to generate example arguments 
## to be stored within serosim and used for examples within the help function.
## This script will follow the same format as README.Rmd
## After the simulation, runserosim outputs will also be stored.
devtools::document("~/Documents/GitHub/serosim")
devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)

set.seed(1)

## Specify the number of time steps in the simulation
times <- seq(1,120,by=1) 
## Set simulation settings argument needed for runserosim 
simulation_settings <- list("t_start"=1,"t_end"=max(times))



## Create an example demography tibble
example_demography <- generate_pop_demography(N=100, times=times, removal_min=0, removal_max=120, prob_removal=0)

## Create example numerical biomarker map for 2 exposures with the same biomarker
example_biomarker_map <- tibble(exposure_id=c("ifxn","vacc"),biomarker_id=c("IgG_titer","IgG_titer"))
## Create an example numerical biomarker map for 2 exposures with the same biomarker
example_biomarker_map_numeric <- tibble(exposure_id=c(1,2),biomarker_id=c(1,1)) 

## Create an example foe_pars array
foe_pars <- array(0, dim=c(1,max(times),n_distinct(example_biomarker_map$exposure_id)))
foe_pars[,,1] <- 0.2
foe_pars[,,2] <- 0.4
example_foe_pars <-foe_pars



## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model (Section 1.5), the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model (Section 1.7) are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
model_pars_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
example_model_pars <- read.csv(file = model_pars_path, header = TRUE)

## Reformat model_pars for runserosim
example_model_pars_numeric<-reformat_biomarker_map(example_model_pars)



## Specify observation_times (serological survey sampling design) to observe biomarker 1 across all individuals at the end of the simulation (t=120)
observation_times<- tibble(i=1:max(example_demography$i),t=120, b=1)


res<- runserosim(
  simulation_settings=simulation_settings,
  demography=example_demography,
  observation_times=observation_times,
  foe_pars=example_foe_pars,
  biomarker_map=example_biomarker_map_numeric,
  model_pars=example_model_pars_numeric,
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
example_exposure_histories_wide<-res$exposure_histories
example_exposure_probabilities<-res$exposure_probabilities_long
example_exposure_force <- res$exposure_force_long
example_biomarker_states<-res$biomarker_states
example_observed_biomarker_states<-res$observed_biomarker_states
example_biomarker_states_wide <- xtabs(value ~ i + t + b, data=example_biomarker_states)

## Additional example data
model_pars_path <- system.file("extdata", "model_pars_typhoid.csv", package = "serosim")
example_model_pars_typhoid <- read.csv(file = model_pars_path, header = TRUE)

## Save all generated data
save(example_demography,file="~/Documents/GitHub/serosim/data/example_demography.rda") 
save(example_biomarker_map,file="~/Documents/GitHub/serosim/data/example_biomarker_map.rda") 
save(example_biomarker_map_numeric,file="~/Documents/GitHub/serosim/data/example_biomarker_map_numeric.rda") 
save(example_foe_pars,file="~/Documents/GitHub/serosim/data/example_foe_pars.rda") 
save(example_model_pars,file="~/Documents/GitHub/serosim/data/example_model_pars.rda") 
save(example_model_pars_numeric,file="~/Documents/GitHub/serosim/data/example_model_pars_numeric.rda") 
save(example_exposure_histories,file="~/Documents/GitHub/serosim/data/example_exposure_histories.rda") 
save(example_exposure_histories_wide,file="~/Documents/GitHub/serosim/data/example_exposure_histories_wide.rda") 
save(example_exposure_probabilities,file="~/Documents/GitHub/serosim/data/example_exposure_probabilities.rda") 
save(example_exposure_force,file="~/Documents/GitHub/serosim/data/example_exposure_force.rda") 
save(example_biomarker_states,file="~/Documents/GitHub/serosim/data/example_biomarker_states.rda") 
save(example_biomarker_states_wide,file="~/Documents/GitHub/serosim/data/example_biomarker_states_wide.rda") 
save(example_observed_biomarker_states,file="~/Documents/GitHub/serosim/data/example_observed_biomarker_states.rda") 
save(example_model_pars_typhoid,file="~/Documents/GitHub/serosim/data/example_model_pars_typhoid.rda") 



