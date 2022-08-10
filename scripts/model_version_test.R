## Create example arguments to test out various model versions

N<-5 
t_start<-1
t_end<-10
times<-t_start:t_end
N_exposure_ids <- 2
N_antigen_ids <-2

## Create demography 
aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),"Group"=list("name"="group","options"=c("1", "2"), "distribution"=c(0.5,0.5)),"NS"=list("name"="NS","options"=c("low", "high"), "distribution"=c(0.5,0.5)) )
demography<-generate_pop_demography(N, times, limit=0, removal_min=2, removal_max=10, prob_removal=0.5, aux=aux)
groups <- demography %>% select(i, group) %>% distinct()
birth_times<-demography %>% select(i, birth) %>% distinct() 
removal_times <- demography %>% select(i, removal) %>% distinct() 

## Create force of infection array 
lambdas<-array(data=1, dim=c(n_distinct(groups$group),max(times),N_exposure_ids))
lambdas[1,1,1]<-3
lambdas[2,2,2]<-2

## Create antigen map 
antigen_map<-tibble(exposure_id=1,antigen_id=1:2)

## Create a tibble with any relevant demographic elements that affect exposure probability 
mod<-tibble(column=c("sex", "sex","NS","NS"), value=c("male","female","low", "high"), modifier=c(1,2,3,4))

## Create a tibble with any relevant age modifiers that affect exposure probability 
age_mod<-tibble(age=0:10, modifier=1:11)

## Establish which part of the loop to work on 
i<-1
t<-1
e<-1
g<-1

i<-2
t<-2
e<-2
g<-2

## Test exposure models
exposure_model_V1(i,t,e, g, lambdas, demography)
exposure_model_V2(i,t,e, g, lambdas, demography, mod=mod)

## The following code gets run earlier within runserosim but it's outputs are needed for the following exposure_model versions
birth_time <- birth_times$birth[i] 
removal_time <- t_end
# simulation_times_tmp <- times[times >= birth_time &  times <= removal_time]

exposure_model_V3(i,t,e, g, lambdas, demography, age_mod)
exposure_model_V4(i,t,e, g, lambdas, demography, mod, age_mod)

## Create dummy arguments to be used for immunity models
N_antigen_ids <-2
antibody_states=array(data=sample(1:10,100,replace=TRUE), dim = c(N,max(times),N_antigen_ids))
vacc_exposures<-c(1,3,5) ## Specify which exposure IDs represent vaccination events 
max_vacc_events<-c(2,1,1,NA,1) ## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
exposure_histories<-array(data=0, dim=c(N,max(times),N_exposure_ids))
exposure_histories[1,2,]<-1

## Import theta
library(readr)
theta <- read.csv("Documents/GitHub/serosim/inst/extdata/theta_V1.csv")

## Test immunity models
immunity_model_V1()
immunity_model_V2(i, t, e, exposure_histories, antibody_states, demography, antigen_map, max_vacc_events)

immunity_model_V3(i, t, e, exposure_histories, antibody_states, demography, antigen_map, theta)
immunity_model_V4(i, t, e, exposure_histories, antibody_states, demography, antigen_map, max_vacc_events, vacc_exposures, theta)


## Test draw_parameters
draw_parameters_V1(i, t, e, demography, theta, antibody_states)
draw_parameters_V2(i, t, e, demography, theta, antibody_states)
draw_parameters_V3(i, t, e, demography, theta, antibody_states)
draw_parameters_V4(i, t, e, demography, theta, antibody_states)

## Create kinetics_parameters to test out antibody models
kinetics_parameters <- vector(mode="list",length=N)
kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                      draw_parameters_V1(i, 2, e, demography, theta, antibody_states))

## Test antibody models
antibody_model_V1(i, t, ag, exposure_histories, kinetics_parameters, antigen_map)

## Reshape antibody_states to test observation model 
antibody_states <- reshape2::melt(antibody_states)
colnames(antibody_states) <- c("i","t","ag","value")
antibody_states <- antibody_states %>% arrange(i, t, ag)

## Create dummy arguments needed for observation models
observation_times <- tibble(i=1:5,t=3, ag=1)
discrete<-c(0,5,8,10) ## Cut offs for discrete assays
boundary<-c(2,10)

## Test out observation models
observation_model_V1(antibody_states, theta, demography, boundary)
observation_model_V2(antibody_states, theta, demography, discrete)
observation_model_V3(antibody_states, theta, demography, boundary)
observation_model_V4(antibody_states, theta, demography, discrete)
  
  
## Reshape runserosim outputs (this is done within runserosim) to test plotting functions
