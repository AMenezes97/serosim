## Create example arguments to test out various model versions

N<-100 
t_start<-1
t_end<-10
times<-t_start:t_end
N_exposure_ids <- 2
N_biomarker_ids <-2

## Create demography 
aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),"Group"=list("name"="group","options"=c("1", "2"), "distribution"=c(0.5,0.5)),"NS"=list("name"="NS","options"=c("low", "high"), "distribution"=c(0.5,0.5)) )
demography<-generate_pop_demography(N, times, limit=0, removal_min=2, removal_max=10, prob_removal=0.5, aux=aux)
groups <- demography %>% select(i, group) %>% distinct()
birth_times<-demography %>% select(i, birth) %>% distinct() 
removal_times <- demography %>% select(i, removal) %>% distinct() 

## Create force of infection array 
foe_pars<-array(data=1, dim=c(n_distinct(groups$group),max(times),N_exposure_ids))
foe_pars[1,1,1]<-3
foe_pars[2,2,2]<-2

## Create biomarker map 
biomarker_map<-tibble(exposure_id=1,biomarker_id=1:2)

## Create a tibble with any relevant demographic elements that affect exposure probability 
dem_mod<-tibble(column=c("sex", "sex","NS","NS"), value=c("male","female","low", "high"), modifier=c(1,2,3,4))

## Create a tibble with any relevant age modifiers that affect exposure probability 
age_mod<-tibble(age=0:10, modifier=1:11)

## Establish which part of the loop to work on 
i<-1
t<-1
x<-1
b<-1

i<-2
t<-2
x<-2
g<-2

## Test exposure models
exposure_model_simple_FOI(i,t,x, g, foe_pars, demography)
exposure_model_dem_mod(i,t,x, g, foe_pars, demography, dem_mod=dem_mod)

## The following code gets run earlier within runserosim but it's outputs are needed for the following exposure_model versions
birth_time <- birth_times$birth[i] 
removal_time <- t_end
# simulation_times_tmp <- times[times >= birth_time &  times <= removal_time]

exposure_model_age_mod(i,t,x, g, foe_pars, demography, age_mod)
exposure_model_dem_age_mod(i,t,x, g, foe_pars, demography, dem_mod, age_mod)

## Create dummy arguments to be used for immunity models
N_biomarker_ids <-1
antibody_states=array(data=sample(1:10,100,replace=TRUE), dim = c(N,max(times),N_biomarker_ids))
vacc_exposures<-c(1,3,5) ## Specify which exposure IDs represent vaccination events 
max_vacc_events<-c(2,1,1,NA,1) ## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
exposure_histories<-array(data=0, dim=c(N,max(times),N_exposure_ids))
exposure_histories[1,2,]<-1

## Import model_pars
library(readr)
## The following file is no longer available
# model_pars <- read.csv("Documents/GitHub/serosim/inst/extdata/model_pars_V1.csv")

## Test immunity models
model_pars_immunity <-tibble(exposure_id=1, biomarker_id=1,name=c("titer_prot_midpoint","titer_prot_width"), mean=c(8,1.5), sd=NA, distribution=NA)
immunity_model_all_successful()
immunity_model_vacc_only(i, t, x, exposure_histories, antibody_states, demography, biomarker_map, max_vacc_events)
immunity_model_vacc_successful_ifxn(i, t, x, exposure_histories, antibody_states, demography, biomarker_map, max_vacc_events, vacc_exposures)
immunity_model_ifxn_titer_prot(i, t, x, exposure_histories, antibody_states, demography, biomarker_map, model_pars=model_pars_immunity)
immunity_model_vacc_ifxn_titer_prot(i, t, x, exposure_histories, antibody_states, demography, biomarker_map, max_vacc_events, vacc_exposures, model_pars=model_pars_immunity)


## Test draw_parameters
draw_parameters_fixed_fx(i, t, x, demography, model_pars, antibody_states)
draw_parameters_random_fx_boost_wane(i, t, x, demography, model_pars, antibody_states)


## Create kinetics_parameters to test out antibody models
## The following lines are run within runserosim
## In order to test the antibody models, we need to draw parameters for exposure
## events prior to the titer time we are trying to calculate within antibody_model
## Create empty kinetics parameters tab 
kinetics_parameters <- vector(mode="list",length=N)
## Draw parameters for an exposure event at time 2
kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                      draw_parameters_fixed_fx(i, 2, x, b, demography, antibody_states, model_pars))

## Test antibody models
antibody_model_biphasic(i, 3, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map)
antibody_model_test(i, 3, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map)


microbenchmark(antibody_model_biphasic(i, 10, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map),
               antibody_model_test(i, 10, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map))


## Reshape antibody_states to test observation model 
antibody_states <- reshape2::melt(antibody_states)
colnames(antibody_states) <- c("i","t","b","value")
antibody_states <- antibody_states %>% arrange(i, t, b)

## Create dummy arguments needed for observation models
observation_times <- tibble(i=1:5,t=3, b=1)
discrete<-c(0,5,8,10) ## Cut offs for discrete assays
boundary<-c(2,10)
model_pars_obs <-tibble(exposure_id=NA, biomarker_id=1:2,name="obs_sd", mean= NA, sd=0.5, distribution="normal")
## Test out observation models
observation_model_continuous_bounded_no_noise(antibody_states, model_pars=model_pars_obs, demography, boundary)
observation_model_discrete_no_noise(antibody_states, model_pars=model_pars_obs, demography, discrete)
observation_model_continuous_bounded_noise(antibody_states, model_pars=model_pars_obs, demography, boundary)
observation_model_discrete_noise(antibody_states, model_pars=model_pars_obs, demography, discrete)

left_join(observation_times,antibody_states)
  
## Reshape runserosim outputs (this is done within runserosim) to test plotting functions
## Antibody states is reshaped above 

## Reshape exposure histories
exposure_histories_long <- NULL
if(sum(exposure_histories) > 0){
  exposure_histories_long <- reshape2::melt(exposure_histories)
  colnames(exposure_histories_long) <- c("i","t","x","value")
  exposure_histories_long <- exposure_histories_long %>% filter(value != 0) %>% select(-value)
  exposure_histories_long <- exposure_histories_long %>% arrange(i, t, x)
}

## Reshape exposure probabilities
exposure_probabilities_long <- reshape2::melt(exposure_probabilities)
colnames(exposure_probabilities_long) <- c("i","t","x","value")
exposure_probabilities_long <- exposure_probabilities_long %>% arrange(i, t, x)

