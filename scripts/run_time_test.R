## This script will run various combinations of models and measure the run time for each 
## Note: This file was before the observation models were updated. Boundary needs to be restructured and changed to bounds
## max_vacc_events might also need to change to max_events


## Set up generic variables needed 
## This script will follow some of the same charcteristic of case study 1:
## 1,000 individuals across 120 time steps
## 2 exposure IDs for 1 antigen ID with 2 observation times 




devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)
library(microbenchmark)

## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))
N<-100

## Pre load the demography categories, values and distributions 
## Specify options for each demography element and the distribution of each within the population
aux <- list("NS"=list("name"="NS","options"=c("low","medium","high"), "distribution"=c(0.3,0.3,0.4)),
            "Group"=list("name"="group","options"=c("1", "2"), "distribution"=c(0.5,0.5)))
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=120, prob_removal=0, aux=aux)

## Set antigen map 
antigen_map <- tibble(exposure_id=c(1,2),antigen_id=c(1,1)) 

## Create an empty array to store the force of infection for all exposure types
foe_pars <- array(0, dim=c(n_distinct(demography$group),max(times),n_distinct(antigen_map$exposure_id)))

## Assign arbitrary forces of infection/vaccination
## Specify the force of vaccination for exposure ID 1 which represents measles natural infection
foe_pars[1,,1] <- 0.2 ## Group 1 (aka Location 1)
foe_pars[2,,1] <- 0.3 ## Group 2 (aka Location 2)

## Specify the force of infection for exposure ID 2 which represents measles vaccination
foe_pars[1,,2] <- 0.4 ## Group 1 (aka Location 1)
foe_pars[2,,2] <- 0.3 ## Group 2 (aka Location 2)


## Specify which exposure IDs represent vaccination events 
vacc_exposures<-2

## Specify the age at which an individual is eligible for vaccination (9 months old for measles)
vacc_age<-c(NA,9)

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
max_vacc_events<-c(NA,1)

## Specify aguements needed for select exposure models
## Create a tibble with any relevant age modifiers that affect exposure probability 
## Individuals who are less than one years old are 3 times more likely be naturally infected
age_mod_1<-tibble(exposure_id=rep(1,11), age=0:10, modifier=c(3,1,1,1,1,1,1,1,1,1,1))

## Individuals who are 0-3 are 2 times more likely to be vaccinated
age_mod_2<-tibble(exposure_id=rep(2,11), age=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))

age_mod<-rbind(age_mod_1,age_mod_2)
age_mod

## Specify the demography exposure modifiers
## Here, individuals who are of low nutritional status are twice as likely of being exposed to diphtheria and pertussis while individuals who are of medium nutritional status are 1.5 times as likely of being exposed when compared to individuals of high nutritional status 
## Individuals of high nutritional status are 3 times more likely to be exposed to exposure ID 1 (vaccination) while individuals who are of medium nutritional status are 2 times more likely to be exposed to exposure ID 1 (vaccination) that individuals of low nutritional status. 
## Note that the modifiers must be defined for all combinations of exposure types and demographic elements
dem_mod<-tibble(exposure_id=c(1,1,1,2,2,2,3,3,3), column=rep("NS",times=9), value=rep(c("low","medium", "high"),3), modifier=c(1,2,3,2,1.5,1,2,1.5,1))
dem_mod

## Specify the number of time steps within a year 
## We are simulating on the monthly scale
t_in_year=12

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model, the titer-ceiling parameters needed for draw_parameters and the observation error needed for the observation model are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
model_pars_path <- system.file("extdata", "model_pars_cs1.csv", package = "serosim")
model_pars <- read.csv(file = model_pars_path, header = TRUE)


## Limits of detection for continuous assays
boundary<-c(20,20000)

## Specify observation_times to observe antigen 1 (aka Measles antibody titer) across all individuals at the midpoint and the end of the simulation (t=60 and t=120)
obs1 <- tibble(i=1:N,t=60, a=1)
obs2 <- tibble(i=1:N,t=120, a=1)
observation_times<-rbind(obs1,obs2)




##Combination 1- (1,1,1,1,1)
Rprof(tmp<-tempfile())

res1<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  antigen_map,
  model_pars,
  exposure_model=exposure_model_simple_FOI,
  immunity_model=immunity_model_all_successful,
  antibody_model=antibody_model_biphasic,
  observation_model=observation_model_continuous_bounded_no_noise,
  draw_parameters=draw_parameters_fixed_fx,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
)

Rprof(NULL)
summaryRprof(tmp)

## Plot antibody states and exposure histories for 10 individuals 
plot_subset_individuals_history(res1$antibody_states, res$exposure_histories_long, subset=10, demography)

## Plot exposure histories for all exposure types
plot_exposure_histories(res1$exposure_histories_long)

## Plot exposure probabilities for all exposure types
plot_exposure_prob(res1$exposure_probabilities_long)

## Plot antibody states for all individuals
plot_titers(res1$antibody_states)

## Plot the first serosurvey at time 60 
obs60<-res1$observed_antibody_states %>% filter(t==60)
plot_obs_titers_one_sample(obs60)

## Plot the second serosurvey at time 120 
obs120<-res1$observed_antibody_states %>% filter(t==120)
plot_obs_titers_one_sample(obs120)

## Plot both serosurveys paired samples
plot_obs_titers_paired_sample(res1$observed_antibody_states)

## Note that the simulated kinetics parameters are also stored
head(res1$kinetics_parameters)

