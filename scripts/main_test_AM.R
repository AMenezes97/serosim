# library(serosim)

##Source the following: 
#source("~/Documents/GitHub/serosim/R/serosim.R")
#source("~/Documents/GitHub/serosim/R/observation_models.R")
#source("~/Documents/GitHub/serosim/R/antibody_models.R")
#source("~/Documents/GitHub/serosim/R/draw_parameters_options.R")
#source("~/Documents/GitHub/serosim/R/immunity_models.R")
#source("~/Documents/GitHub/serosim/R/exposure_models.R")
#source("~/Documents/GitHub/serosim/R/generate_pop_demography.R")
#source("~/Documents/GitHub/serosim/R/generate_plots.R")
  
library(tidyverse)
library(data.table)
library(ggplot2)

##******************Component 1: Simulation Settings**************************** 
## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))



##******************Component 2: Population Demography************************** 
## Specify the number of individuals in the simulation 
N <- 100

## Pre load the demography categories, values and distributions 
aux <- list("SES"=list("name"="SES","options"=c("low","medium","high"), "distribution"=c(0.2,0.2,0.6)),
            "NS"=list("name"="NS","options"=c("low","high"),"distribution"=c(0.5,0.5)),
            "Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),
            "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "distribution"=c(0.25,0.25,0.25,0.25)) )

## Simulate demography settings 
demography <- generate_pop_demography(N, times, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3, aux=aux)



##*****************Component 3: Force of Infection and Exposure Model***********
## Specify the number of exposures (eventually, this won't be needed because the serosim function can extract this from the biomarker_map)
N_exposure_ids <- 2

## Specify the force of infection array (different format)
group<-c(1,2,3,4)
foe_pars <- array(rep(0.5,length(group)), dim=c(length(group),max(times),N_exposure_ids))

## Specify exposure model within serosim function below

## Some exposure models will require the following arguments specified within serosim

## Create a tibble with any relevant demographic elements that affect exposure probability 
# mod<-tibble(column=c("sex", "sex","NS","NS"), value=c("male","female","low", "high"), modifier=c(1,2,3,4))

## Create a tibble with any relevant age modifiers that affect exposure probability 
# age_mod<-tibble(age=0:10, modifier=1:11)



##***********************Component 4: Biomarker Map*******************************
## Specify which biomarkers are present in each exposure type
biomarker_map <- tibble(exposure_id=c(1,1,2),biomarker_id=c(1,2,1)) 



##**********************Component 5: Immunity Model*****************************
## Specify immunity model within serosim function below 

## Some immunity models will require the following arguments specified within serosim:

## Specify which exposure IDs represent vaccination events 
vacc_exposures<-1

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
max_vacc_events<-1 




##****Component 6: Antibody Model, Antibody Kinetics Parameters, and draw_parameters*****
## Specify antibody model within serosim function below 

## Specify antibody kinetics parameters 
model_pars <- read.csv("inst/extdata/model_pars_test_1.csv")

## Specify draw_parameters within serosim function below 




##*************Component 7: Observation Model and observation_times*************
## Specify observation model within serosim function below 

## Some observation models will require the following arguments specified within serosim:

## Cut offs for discrete assays
# discrete<-c(0,5,8,10) 

## Limits of detection for continuous assays
boundary<-c(2,20)

## Set observation settings 
obs1 <- tibble(i=1:N,t=60, b=1)
obs2 <- tibble(i=1:N,t=60, b=2)
obs3 <- tibble(i=1:N,t=120, b=1)
obs4 <- tibble(i=1:N,t=120, b=2)
observation_times<-rbind(obs1,obs2,obs3,obs4)

    
##***************************Run Simulation*************************************
# Rprof(tmp<-tempfile())
## Test full function with generated inputs
res<- runserosim(
  simulation_settings,
  demography, 
  observation_times,
  foe_pars, 
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_simple_FOE, 
  immunity_model=immunity_model_vacc_ifxn_titer_prot, 
  antibody_model=antibody_model_biphasic, 
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_fixed_fx, 
  
  ## Pre-specified parameters/events
  exposure_histories_fixed=NULL,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures
)

# Rprof(NULL)
# summaryRprof(tmp)


# # ## Generate Plots 
# ggplot(res$antibody_states) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~b)
plot_titers(res$antibody_states)
# ggplot(res$exposure_probabilities_long) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~e)
plot_exposure_prob(res$exposure_probabilities_long)
# ggplot(res$observed_antibody_states) + geom_jitter(aes(x=t,y=observed),height=0.1,width=0.25) + facet_wrap(~b) + scale_x_continuous(limits=range(times))
plot_obs_titers_one_sample(res$observed_antibody_states)
plot_obs_titers_paired_sample(res$observed_antibody_states)
plot_exposure_histories(res$exposure_histories_long)



