devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)

## Case Study 2
times <- seq(1,120,by=1) 
simulation_settings <- list("t_start"=1,"t_end"=max(times))
N<-10
aux <- list("NS"=list("name"="NS","options"=c("low","medium","high"), "distribution"=c(0.3,0.3,0.4)),
            "Group"=list("name"="group","options"=c("1", "2"), "distribution"=c(0.5,0.5)))
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=120, prob_removal=0.2, aux=aux)
biomarker_map_original <- tibble(exposure_id=c("AB_vacc","AB_vacc","A_ifxn","B_ifxn"),biomarker_id=c("A_antibody","B_antibody","A_antibody","B_antibody"))
biomarker_map <-reformat_biomarker_map(biomarker_map_original)
foe_pars <- array(0, dim=c(n_distinct(demography$group),max(times), n_distinct(biomarker_map$exposure_id)))
foe_pars[1,,1] <- 0.2 ## Group 1 (aka Group 1)
foe_pars[2,,1] <- 0.3 ## Group 2 (aka Group 2)
foe_pars[1,,2] <- 0.4 ## Group 1 (aka Group 1)
foe_pars[2,,2] <- 0.3 ## Group 2 (aka Group 2)
foe_pars[1,,3] <- 0.2 ## Group 1 (aka Group 1)
foe_pars[2,,3] <- 0.1 ## Group 2 (aka Group 2)
age_mod_1<-tibble(exposure_id=rep(1,11), age=0:10, modifier=c(3,1,1,1,1,1,1,1,1,1,1))
age_mod_2<-tibble(exposure_id=rep(2,11), age=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))
age_mod_3<-tibble(exposure_id=rep(3,11), age=0:10, modifier=c(1,2,1,1,1,1,1,1,1,1,1))
age_mod<-rbind(age_mod_1,age_mod_2,age_mod_3)
dem_mod<-tibble(exposure_id=c(1,1,1,2,2,2,3,3,3), column=rep("NS",times=9), value=rep(c("low","medium", "high"),3), modifier=c(1,2,3,2,1.5,1,2,1.5,1))
t_in_year=12
vacc_exposures<-1
vacc_age<-c(2,NA,NA)
max_vacc_events<-c(1,NA,NA)
model_pars_path <- system.file("extdata", "model_pars_cs2.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars<-reformat_model_pars(biomarker_map_original,model_pars_original )
bounds<-tibble(biomarker_id=c(1,1,2,2),name=rep(c("lower_bound","upper_bound"),2),value=c(2,50,4,100))
obs1 <- tibble(i=1:N,t=120, b=1)
obs2 <- tibble(i=1:N,t=120, b=2)
observation_times<-rbind(obs1,obs2)
cutoffs<-matrix(data=0,nrow = 2,ncol = 5)
biomarker1<-c(0,2,5,10,15)
biomarker2<-c(0,20,50,75,100)
cutoffs[1,]<-biomarker1
cutoffs[2,]<-biomarker2

res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_age_mod,
  immunity_model=immunity_model_vacc_ifxn_titer_prot,
  antibody_model=antibody_model_biphasic,
  observation_model=observation_model_continuous,
  draw_parameters=draw_parameters_random_fx_titer_dep,

  ## Other arguments needed
  bounds=bounds,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  dem_mod=dem_mod,
  age_mod=age_mod,
  t_in_year=t_in_year,
  cutoffs=cutoffs
)




## Plot antibody states and exposure histories for 10 individuals 
plot_subset_individuals_history(res$antibody_states, res$exposure_histories_long, subset=10, demography)
## Plot exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)
## Plot exposure probabilities for all exposure types
plot_exposure_prob(res$exposure_probabilities_long)
## Plot antibody states for all individuals
plot_titers(res$antibody_states)
## Plot the serosurvey results 
plot_obs_titers_one_sample(res$observed_antibody_states)

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)
