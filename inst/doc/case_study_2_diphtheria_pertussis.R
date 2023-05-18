## ---- echo=FALSE--------------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## ----message=FALSE, warning=FALSE,eval=TRUE-----------------------------------
## Install and load serosim 
## devtools::install_github("AMenezes97/serosim")
library(serosim)

## Load additional packages required 
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)

## -----------------------------------------------------------------------------
## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))

## -----------------------------------------------------------------------------
## Specify the number of individuals in the simulation 
N<-100

## Pre-load the demography categories, values and distributions 
## Specify options for each demography element and the distribution of each within the population
## We are interested in simulating a population where individuals have varying nutritional status and can reside in either of 2 locations
aux <- list("NS"=list("name"="NS","options"=c("low","medium","high"), "distribution"=c(0.3,0.3,0.4)),
            "Group"=list("name"="group","options"=c(1, 2), "distribution"=c(0.5,0.5)))


## Generate the population demography tibble
## Let's assume that individuals are removed from the population and set prob_removal to 0.2
demography <- generate_pop_demography(N, times, age_min=0, removal_min=1, removal_max=120, prob_removal=0.2, aux=aux)

## Examine the generated demography tibble
head(demography)
tail(demography)

## -----------------------------------------------------------------------------
## Create biomarker map
biomarker_map_original <- tibble(exposure_id=c("DP_ifxn","PT_ifxn","vacc","vacc"),biomarker_id=c("DP_antibody","PT_antibody","DP_antibody","PT_antibody"))
biomarker_map_original

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map_original)
biomarker_map

## ---- fig.dim = c(5, 6)-------------------------------------------------------
## Create an empty array to store the force of exposure for all exposure types across all time steps and groups
## Dimension 1: Group
## Dimension 2: Times
## Dimension 3: Exposure ID in the numeric order that they appear in the biomarker map
foe_pars <- array(0, dim=c(n_distinct(demography$group),max(times), n_distinct(biomarker_map$exposure_id)))

## Note that we can specify a different force of exposure for each group, time and exposure ID
## We specified the same value for all time steps within foe_pars for simplicity.

## Specify the force of exposure for exposure ID 1 which represents diphtheria natural infection (DP_ifxn)
foe_pars[1,,1] <- 0.04 ## Group 1 (aka Location 1)
foe_pars[2,,1] <- 0.03 ## Group 2 (aka Location 2)

## Specify the force of exposure for exposure ID 2 which represents pertussis natural infection (PT_ifxn)
foe_pars[1,,2] <- 0.02 ## Group 1 (aka Location 1)
foe_pars[2,,2] <- 0.01 ## Group 2 (aka Location 2)

## Specify the force of exposure for exposure ID 3 which represents diphtheria and pertussis combined vaccine (vacc)
foe_pars[1,,3] <- 0.02 ## Group 1 (aka Location 1)
foe_pars[2,,3] <- 0.03 ## Group 2 (aka Location 2)

## Specify a simple exposure model which calculates the probability of exposure from the force of exposure modulated by age and demography elements 
exposure_model<-exposure_model_dem_mod

## This exposure model requires dem_mod and t_in_year arguments

## Create a tibble with any relevant demography modifiers that affect exposure probability 
## For simplicity, we selected arbitrary numbers. 

## First, we will specify age modifiers 
## Individuals who are 0-3 are 2 times more likely to be exposed to diphtheria
age_mod_1<-tibble(exposure_id=rep(1,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))

## Individuals who are 0-3 are 2 times more likely to be exposed to pertussis
age_mod_2<-tibble(exposure_id=rep(2,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,1,1,1,1,1,1,1,1))

## Individuals who are less than one year old are 3 times more likely be vaccinated than other age classes
age_mod_3<-tibble(exposure_id=rep(3,11), column=rep("age",times=11),  value=0:10, modifier=c(3,1,1,1,1,1,1,1,1,1,1))
age_mod<-rbind(age_mod_1,age_mod_2,age_mod_3)
age_mod

## Now we will  specify additional demography exposure modifiers and combine them with the previous ones
## Here, individuals who are of low nutritional status are twice as likely of being exposed to diphtheria and pertussis while individuals who are of medium nutritional status are 1.5 times as likely of being exposed when compared to individuals of high nutritional status 
## Individuals of high nutritional status are 3 times more likely to be exposed to exposure ID 3 (vaccination) while individuals who are of medium nutritional status are 2 times more likely to be exposed to exposure ID 3 (vaccination) than individuals of low nutritional status. 
## Note that the modifiers must be defined for all combinations of exposure types and demographic elements
mod<-tibble(exposure_id=c(1,1,1,2,2,2,3,3,3), column=rep("NS",times=9), value=rep(c("low","medium", "high"),3), modifier=c(2,1.5,1,2,1.5,1,1,2,3))
mod
## Combine both age modifiers and additional modifiers
dem_mod<-rbind(age_mod,mod)

## Specify the number of time steps within a year which will be used to calculate an individual's age. 
## We are simulating on the monthly scale so there are 12 time steps within a year.
t_in_year=12

## Examine the probability of exposure over time for the specified exposure model
plot_exposure_model(indivs=1:5,exposure_model=exposure_model_dem_mod, times=times, n_groups = 2, n_exposures = 3, foe_pars=foe_pars, demography=demography, dem_mod=dem_mod,t_in_year=t_in_year)

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Specify immunity model within the runserosim function below 
immunity_model<-immunity_model_vacc_ifxn_biomarker_prot

## Specify which exposure IDs represent vaccination events 
## In this case, only exposure ID 3 is a vaccination event
vacc_exposures<-3

## Specify the time step after birth at which an individual is eligible for vaccination (2 months old for diphtheria and pertussis combine vaccine); ; note non-vaccine exposures are listed as NAs
vacc_age<-c(NA,NA,2)

## Specify the maximum number of successful exposure events an individual can receive for each exposure type
## We placed no successful exposure limit on the number of diphtheria or pertussis infection exposures and 3 dose limit on the vaccine exposure.
max_events<-c(Inf,Inf,3)

## Plot biomarker-mediated protection curve given parameters specified within model_pars for biomarker 1, diphtheria antibody  (DP_antibody) which will be loaded in section 1.6
## The immunity model we selected will take into account an individual's current titer to diphtheria when determining the probability of successful infection. 
## Titers are ploted in mIU/mL 
plot_biomarker_mediated_protection(0:150, biomarker_prot_midpoint=75, biomarker_prot_width=.1)

## Plot biomarker-mediated protection curve given parameters specified within model_pars for biomarker 2, pertussis antibody (PT_antibody) which will be loaded in section 1.6
## The immunity model we selected will take into account an individual's current titer to pertussis when determining the probability of successful infection. 
## Titers are ploted in IU/mL 
plot_biomarker_mediated_protection(0:100, biomarker_prot_midpoint=40, biomarker_prot_width=.25)

## -----------------------------------------------------------------------------
## Specify the antibody model 
antibody_model<-antibody_model_biphasic

## Bring in the antibody parameters needed for the antibody model 
## Note that the titer-mediated protection parameters needed for the immunity model, the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model are all defined here too.
## Also note that these are all arbitrary parameter values chosen for this toy example.
model_pars_path <- system.file("extdata", "model_pars_cs2.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars_original 

## Reformat model_pars for runserosim
model_pars<-reformat_biomarker_map(model_pars_original)
model_pars

## Specify the draw_parameters function to use 
draw_parameters<-draw_parameters_random_fx_biomarker_dep

## ---- fig.dim = c(5, 4)-------------------------------------------------------
## Plot biomarker (in this case antibody titer) dependent boosting effects given parameters specified within model_pars for biomarker 1, diphtheria antibody (DP_antibody)
plot_biomarker_dependent_boosting(start=0, end=2, by=.1, biomarker_ceiling_threshold=1.7, biomarker_ceiling_gradient=0.529411)

## Plot biomarker (in this case antibody titer) dependent boosting effects given parameters specified within model_pars for biomarker 2, pertussis antibody (PT_antibody)
plot_biomarker_dependent_boosting(start=0, end=125, by=1, biomarker_ceiling_threshold=100, biomarker_ceiling_gradient=0.009)

## Plot example biomarker trajectories given the specified antibody kinetics model, model parameters and draw parameters function 
plot_antibody_model(antibody_model_biphasic, N=25, model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx_biomarker_dep, biomarker_map=biomarker_map)

## -----------------------------------------------------------------------------
## Specify the limits of detection for each biomarker for the continuous assays
bounds<-tibble(biomarker_id=c(1,1,2,2),name=rep(c("lower_bound","upper_bound"),2),value=c(0.01,2,5,200))

## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify observation_times to observe both biomarkers (aka DP_antibody and PT_antibody titers) across all individuals at the end of the simulation (t=120)
obs1 <- tibble(i=1:N,t=120, b=1)
obs2 <- tibble(i=1:N,t=120, b=2)
observation_times<-rbind(obs1,obs2)

## ----message=FALSE, warning=FALSE, eval=TRUE----------------------------------
res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model,
  immunity_model,
  antibody_model,
  observation_model,
  draw_parameters,

  ## Other arguments needed
  bounds=bounds,
  max_events=max_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  dem_mod=dem_mod,
  t_in_year=t_in_year,
  VERBOSE=NULL,
  attempt_precomputation = FALSE
)
## Note that models and arguments specified earlier in the code can be specified directly within this function.

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot diphtheria biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states %>% filter(b==1), res$exposure_histories_long %>% filter(x==1 | x==3), subset=10, demography, removal=TRUE)

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot pertussis biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states %>% filter(b==2), res$exposure_histories_long %>% filter(x==2 | x==3), subset=10, demography, removal=TRUE)

## ---- fig.dim = c(5, 6)-------------------------------------------------------
## Plot individual probability of exposure for all exposure types.
## This is the output of the exposure model.
plot_exposure_force(res$exposure_force_long)

## ---- fig.dim = c(5, 6)-------------------------------------------------------
## Plot individual successful exposure probabilities for all exposure types
## This is the output of the exposure model multiplied by the output of the immunity model.
## In other words, this is the probability of exposure event being successful and inducing an immunological response
plot_exposure_prob(res$exposure_probabilities_long)

## ---- fig.dim = c(5, 6)-------------------------------------------------------
## Plot individual exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot diphtheria antibody (biomarker 1) states for all individuals (true biomarker quantities)
plot_biomarker_quantity(res$biomarker_states %>% filter(b==1))

## ---- , fig.dim = c(4, 4.5)---------------------------------------------------
## Plot pertussis antibody (biomarker 2) states for all individuals (true biomarker quantities)
plot_biomarker_quantity(res$biomarker_states %>% filter(b==2))

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot the diphtheria serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states %>% filter(b==1))

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot the pertussis serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states %>% filter(b==2))

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)

