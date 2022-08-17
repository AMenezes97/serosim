# library(serosim)

##Source the following: 
## serosim.R
## generate_pop_demography.R
## exposure_models.R
## immunity_models.R
## draw_parameters_versions.R
## antibody_models.R
## oobservation_models.R
  


##******************Component 1: Simulation Settings**************************** 
## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))



##******************Component 2: Population Demography************************** 
## Specify the number of individuals in the simulation 
N <- 500

## Pre load the demography categories, values and distributions 
aux <- list("SES"=list("name"="SES","options"=c("low","medium","high"), "distribution"=c(0.2,0.2,0.6)),
            "NS"=list("name"="NS","options"=c("low","high"),"distribution"=c(0.5,0.5)),
            "Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),
            "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "distribution"=c(0.25,0.25,0.25,0.25)) )

## Simulate demography settings 
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=120, prob_removal=0.3, aux=aux)



##*****************Component 3: Force of Infection and Exposure Model***********
## Specify the number of exposures (eventually, this won't be needed because the serosim function can extract this from the antigen_map)
N_exposure_ids <- 1

## Specify the force of infection array (different format)
group<-c(1,2,3,4)
lambdas <- array(rep(0.5,length(group)), dim=c(length(group),max(times),N_exposure_ids))

## Specify exposure model within serosim function below

## Some exposure models will require the following arguments specified within serosim

## Create a tibble with any relevant demographic elements that affect exposure probability 
# mod<-tibble(column=c("sex", "sex","NS","NS"), value=c("male","female","low", "high"), modifier=c(1,2,3,4))

## Create a tibble with any relevant age modifiers that affect exposure probability 
# age_mod<-tibble(age=0:10, modifier=1:11)



##***********************Component 4: Antigen Map*******************************
## Specify which antigens are present in each exposure type
antigen_map <- tibble(exposure_id=1,antigen_id=c(1,2)) 



##**********************Component 5: Immunity Model*****************************
## Specify immunity model within serosim function below 

## Some immunity models will require the following arguments specified within serosim:

## Specify which exposure IDs represent vaccination events 
# vacc_exposures<-c(1,3,5) 

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
# max_vacc_events<-c(2,1,1,NA,1) 




##****Component 6: Antibody Model, Antibody Kinetics Parameters, and draw_parameters*****
## Specify antibody model within serosim function below 

## Specify antibody kinetics parameters 
theta <- read.csv("Documents/GitHub/serosim/inst/extdata/theta_V1.csv")

## Specify draw_parameters within serosim function below 




##*************Component 7: Observation Model and observation_times*************
## Specify observation model within serosim function below 

## Some observation models will require the following arguments specified within serosim:

## Cut offs for discrete assays
# discrete<-c(0,5,8,10) 

## Limits of detection for continuous assays
boundary<-c(2,10)

## Set observation settings 
observation_times <- tibble(i=1:N,t=120, ag=1)



    
##***************************Run Simulation*************************************

## Test full function with generated inputs
res<- serosim(
  simulation_settings,
  demography, 
  observation_times,
  lambdas, 
  antigen_map,
  theta,
  exposure_model=exposure_model_simple_FOI, 
  immunity_model=immunity_model_all_successful, 
  antibody_model=antibody_model_biphasic, 
  observation_model=observation_model_continuous_bounded_no_noise,
  draw_parameters=draw_parameters_fixed_fx, 
  
  ## Pre-specified parameters/events
  exposure_histories_fixed=NULL,
  
  ## Other arguments needed
  boundary=boundary
)


## Generate Plots 
ggplot(res$antibody_states) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~ag)
ggplot(res$exposure_probabilities_long) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~e)
ggplot(res$observed_antibody_states) + geom_jitter(aes(x=t,y=value),height=0.1,width=0.25) + facet_wrap(~ag) + scale_x_continuous(limits=range(times))


## Examine serosim outputs 
res$exposure_histories
res$exposure_histories_long ## Create plot 
res$exposure_probabilities
res$exposure_probabilities_long ## Create plot
res$antibody_states
res$observed_antibody_states ## Create plot 
res$kinetics_parameters

