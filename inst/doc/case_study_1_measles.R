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

## Generate the population demography tibble
## See help file(?generate_pop_demography) for more information on function arguments.
## age_min is set to 0 which allows births to occur until the last time step
## Let's assume that no individuals are removed from the population and set prob_removal to 0
demography <- generate_pop_demography(N, times, age_min=0, prob_removal=0)

## Examine the generated demography tibble
summary(demography)

## -----------------------------------------------------------------------------
## Create biomarker map
biomarker_map_original <- tibble(exposure_id=c("measles_ifxn","measles_vacc"),biomarker_id=c("measles_IgG","measles_IgG"))
biomarker_map_original

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map_original)
biomarker_map

## ---- fig.dim = c(4, 3)-------------------------------------------------------
## Create an empty array to store the force of exposure for all exposure types
## Dimension 1: Group
## Dimension 2: Times
## Dimension 3: Exposure ID 
foe_pars <- array(0, dim=c(1,max(times),n_distinct(biomarker_map$exposure_id)))

## Specify the force of exposure for exposure ID 1 which represents measles natural infection
foe_pars[,,1] <- 0.02

## Specify the force of exposure for exposure ID 2 which represents measles vaccination
foe_pars[,,2] <- 0.04

## Specify a simple exposure model which calculates the probability of exposure directly from the force of exposure at that time step
exposure_model<-exposure_model_simple_FOE
## In this selected model, the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.

## Examine the probability of exposure over time for the specified exposure model
plot_exposure_model(exposure_mode=exposure_model_simple_FOE, times=times, n_groups = 1, n_exposures = 2, foe_pars=foe_pars)

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Specify immunity model within runserosim function below 
immunity_model<-immunity_model_vacc_ifxn_biomarker_prot
## This immunity model requires 3 additional arguments specified below. 

## Specify which exposure IDs represent vaccination events 
vacc_exposures<-2

## Specify the age at which an individual is eligible for measles vaccination (9 months old); note non vaccine exposures are listed as NAs
vacc_age<-c(NA,9)

## Specify the maximum number of successful exposure events an individual can receive for each exposure type
## For this example, we are assuming that there is only one dose of the measles vaccine and one possible measles infection
max_events<-c(1,1)

## Plot biomarker-mediated protection curve (also known as titer-mediated protection) given parameters specified within model_pars for biomarker 1 (measles_IgG) which will be loaded in section 1.6
plot_biomarker_mediated_protection(0:1000, biomarker_prot_midpoint=200, biomarker_prot_width=1)

## -----------------------------------------------------------------------------
## Specify the antibody model 
antibody_model<-antibody_model_biphasic

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model (Section 1.5), the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model (Section 1.7) are all defined here too.
model_pars_path <- system.file("extdata", "model_pars_cs1.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars_original 

## Reformat model_pars for runserosim
model_pars<-reformat_biomarker_map(model_pars_original)
model_pars

## Specify the draw_parameters function to use 
draw_parameters<-draw_parameters_random_fx_biomarker_dep

## Plot biomarker dependent boosting effects given parameters specified within model_pars for biomarker 1 (measles_IgG)
plot_biomarker_dependent_boosting(start=0, end=1000, by=1, biomarker_ceiling_threshold=200, biomarker_ceiling_gradient=.0045)

## These parameters indicate that any individual with a titer level above 200 mIU/mL at the time of a second exposure event will only receive 10% of the regular boosting parameter for that event. Individuals with a measles titer of 0 will receive the full boost.

## ---- fig.dim = c(5, 4)-------------------------------------------------------
## Plot example biomarker trajectories given the specified antibody kinetics model, model parameters and draw parameters function 
plot_antibody_model(antibody_model_biphasic, N=100, model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx_biomarker_dep, biomarker_map=biomarker_map)

## -----------------------------------------------------------------------------
## Specify the limits of detection for each biomarker for the continuous assays
bounds<-tibble(biomarker_id=c(1,1),name=c("lower_bound","upper_bound"),value=c(100,Inf))

## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify assay sensitivity and specificity needed for the observation model
sensitivity<-0.996
specificity<-1

## Specify observation_times (serological survey sampling design) to observe biomarker 1 (measles IgG antibody titer) across all individuals at the midpoint and the end of the simulation (t=60 and t=120)
obs1 <- tibble(i=1:N,t=60, b=1)
obs2 <- tibble(i=1:N,t=120, b=1)
observation_times<-rbind(obs1,obs2)

## -----------------------------------------------------------------------------
## If VERBOSE is specified within runserosim, there will be an update message printed for every multiple of the integer specified.
## Print a message for every 10 individuals 
VERBOSE<-10

## ----message=FALSE, warning=FALSE, eval=TRUE----------------------------------
## Run the core simulation and save outputs in "res"
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

  ## Specify other arguments needed
  VERBOSE=VERBOSE,
  bounds=bounds,
  max_events=max_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  sensitivity=sensitivity,
  specificity=specificity
)

## Note that models and arguments specified earlier in the code can be specified directly within this function.

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states, res$exposure_histories_long, subset=10, demography)

## ---- fig.dim = c(4, 6)-------------------------------------------------------
## Plot individual probability of exposure for all exposure types.
## This is the output of the exposure model.
## Note: All individuals are under the same force of exposure since we specified a simple exposure model and constant foe_pars
plot_exposure_force(res$exposure_force_long)

## ----fig.dim = c(4, 6)--------------------------------------------------------
## Plot individual successful exposure probabilities for all exposure types
## This is the output of the exposure model multiplied by the output of the immunity model.
## In other words, this is the probability of exposure event being successful and inducing an immunological response
plot_exposure_prob(res$exposure_probabilities_long)

## ----fig.dim = c(4, 6)--------------------------------------------------------
## Plot individual exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)

## ---- fig.dim = c(4, 6)-------------------------------------------------------
## Plot true biomarker quantities for all individuals across the entire simulation period
plot_biomarker_quantity(res$biomarker_states)

## ----fig.dim = c(4, 5)--------------------------------------------------------
## Plot the serosurvey results (observed biomarker quantities) at time 60 
obs60<-res$observed_biomarker_states %>% filter(t==60)
plot_obs_biomarkers_one_sample(obs60)

## ---- fig.dim = c(4, 5)-------------------------------------------------------
## Plot the serosurvey results (observed biomarker quantities) at time 120 
obs120<-res$observed_biomarker_states %>% filter(t==120)
plot_obs_biomarkers_one_sample(obs120)

## ----fig.dim = c(5, 5)--------------------------------------------------------
## Plot both serosurveys paired samples
plot_obs_biomarkers_paired_sample(res$observed_biomarker_states)

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)


