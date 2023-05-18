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
biomarker_map_original <- tibble(exposure_id=c("strain_a","strain_a","strain_a","strain_b","strain_b","strain_b", "strain_c","strain_c","strain_c"),biomarker_id=rep(c("biomarker_a","biomarker_b","biomarker_c"),3))
biomarker_map_original

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map_original)
biomarker_map

## -----------------------------------------------------------------------------
## Create an empty array to store the force of exposure for all exposure types
foe_pars <- array(0, dim=c(1,max(times),n_distinct(biomarker_map$exposure_id)))
## Specify the force of exposure for exposure ID 1 which represents strain A circulation
foe_pars[,1:40,1] <- 0.1
foe_pars[,41:80,1] <- 0
foe_pars[,81:120,1] <- 0
## Specify the force of exposure for exposure ID 2 which represents strain B circulation
foe_pars[,1:40,2] <- 0
foe_pars[,41:80,2] <- 0.1
foe_pars[,81:120,2] <- 0
## Specify the force of exposure for exposure ID 3 which represents strain C circulation
foe_pars[,1:40,3] <- 0
foe_pars[,41:80,3] <- 0
foe_pars[,81:120,3] <- 0.1

## Specify a simple exposure model which calculates the probability of exposure directly from the force of exposure at that time step. In this selected model, the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.
exposure_model<-exposure_model_simple_FOE

## ---- fig.dim = c(6, 3)-------------------------------------------------------
## Examine the probability of exposure to each strain over time for the specified exposure model and foe_pars array
plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times, n_groups = 1, n_exposures = 3, foe_pars=foe_pars)

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Specify immunity model within the runserosim function below 
immunity_model<-immunity_model_ifxn_biomarker_prot
## The immunity model we selected will take into account an individual's current biomarker quantity to all three strains and the cross-reactivity between different biomarkers when determining the probability of successful infection. 

## Specify the maximum number of successful exposure events an individual can receive for each exposure type
## We placed no successful exposure limit on the number influenza infection exposures
max_events<-c(Inf,Inf,Inf)

## Specify the cross-reactivity table argument which will be used by the immunity model. 
cr_table <-tibble(
  exposure_id=c("strain_a","strain_a","strain_a","strain_b","strain_b","strain_b", "strain_c","strain_c","strain_c"),
  biomarker_id=rep(c("biomarker_a","biomarker_b","biomarker_c"),3),
  cross_reactivity=c(1,0.5,0.25,0.5,1,0.5,0.25,.5,1))

## The cross_reactivity table also requires that exposure_id and biomarker_id are numeric. The reformat_biomarker_map function will convert both columns to numbers. 
cr_table<-reformat_biomarker_map(cr_table)

## ---- fig.dim = c(3, 3)-------------------------------------------------------
## Plot biomarker-mediated protection curve given parameters specified within model_pars for biomarker 1, influenza strain A antibody  (DP_antibody) which will be loaded in section 1.6. Biomarker mediated protection parameters are the same for all 3 influenza strains so the following plot will look the same for all exposure/biomarker types.
plot_biomarker_mediated_protection(0:120, biomarker_prot_midpoint=40, biomarker_prot_width=.09) 

## -----------------------------------------------------------------------------
## Specify the antibody model 
antibody_model<-antibody_model_monophasic

## Bring in the antibody parameters needed for the antibody model
## Note that the observation error parameter needed for the observation model (Section 1.7) is defined here too.
model_pars_path <- system.file("extdata", "model_pars_cs3.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars_original

# model_pars_original <- read.csv("~/Documents/GitHub/serosim/inst/extdata/model_pars_cs3.csv")
# model_pars_original

## Reformat model_pars for runserosim
model_pars<-reformat_biomarker_map(model_pars_original)
model_pars

## Specify the draw_parameters function
draw_parameters<-draw_parameters_random_fx

## ---- fig.dim = c(8, 4)-------------------------------------------------------
## Plot example biomarker trajectories given the specified antibody kinetics model, model parameters and draw parameters function 
plot_antibody_model(antibody_model_monophasic, N=100, model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx, biomarker_map=biomarker_map)

## -----------------------------------------------------------------------------
## Specify the observation model 
observation_model<-observation_model_discrete

## Specify the cutoffs for the discrete assay. This will depend on the dilutions of the haemagglutination assay.
## This is a matrix with each row containing all of the cutoffs for that biomarker. Here, we have set the same cutoffs for all three biomarkers 
breaks <- seq(0,120,by=10)
cutoffs <- matrix(breaks,nrow=3,ncol=length(breaks), byrow=TRUE)

## Specify assay sensitivity and specificity needed for the observation model
sensitivity<-0.85
specificity<-0.9

## Specify observation_times (serological survey sampling design) to observe all three biomarkers across all individuals at the end of the simulation (t=120)
obs1 <- tibble(i=1:N,t=120, b=1)
obs2 <- tibble(i=1:N,t=120, b=2)
obs3 <- tibble(i=1:N,t=120, b=3)
observation_times<-rbind(obs1,obs2,obs3)

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

  ## Other arguments needed
  max_events=max_events,
  cutoffs=cutoffs,
  cross_reactivity_table=cr_table,
  sensitivity=sensitivity,
  specificity=specificity
)

## Note that models and arguments specified earlier in the code can be specified directly within this function.

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot influenza strain A specific biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states %>% filter(b==1), res$exposure_histories_long, subset=10, demography, removal=TRUE)

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot influenza strain B specific biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states %>% filter(b==2), res$exposure_histories_long, subset=10, demography, removal=TRUE)

## ---- fig.dim = c(6, 4)-------------------------------------------------------
## Plot influenza strain C specific biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states %>% filter(b==3), res$exposure_histories_long, subset=10, demography, removal=TRUE)

## ---- fig.dim = c(5, 6)-------------------------------------------------------
## Plot individual probability of exposure for all exposure types.
## This is the output of the exposure model. Note that this reflects our inputs that the influenza strains were never co-circulating.
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
## Plot influenza strain A antibody (biomarker 1) states for all individuals (true biomarker quantities)
plot_biomarker_quantity(res$biomarker_states %>% filter(b==1))

## ---- , fig.dim = c(4, 4.5)---------------------------------------------------
## Plot influenza strain B antibody (biomarker 2) states for all individuals (true biomarker quantities)
plot_biomarker_quantity(res$biomarker_states %>% filter(b==2))

## ---- , fig.dim = c(4, 4.5)---------------------------------------------------
## Plot influenza strain C antibody (biomarker 3) states for all individuals (true biomarker quantities)
plot_biomarker_quantity(res$biomarker_states %>% filter(b==3))

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot the influenza strain A serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states %>% filter(b==1))

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot the influenza strain B serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states %>% filter(b==2))

## ---- fig.dim = c(4, 4.5)-----------------------------------------------------
## Plot the influenza strain C serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states %>% filter(b==3))

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)

