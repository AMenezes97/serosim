## Simulate a data set to input to Keya's latent class model 

## Data structure
## mom_indicator:  binary indicator of if an individual had mother recall of any mcv dose among children with or without card
## card_indicator: binary indicator of if an individual had record (via card) of any mcv dose
## none_indicator = binary indicator of if an individual had neither card nor recall
## boost_indicator = binary indicator of if an individual has a suspected boost
## titer = log antibody concentrations of each individual 
## age = age in years of each individual. Age here is a double (e.g. 4.39) (Age ranges 0.7 to 14.99)
## ageyr = age in years of each individual as an integer (floor of age) 
## ageyr_obs = age in integers among those with observed vaccine information
## ageyr_mis = age in integers among those with missing vaccine information 
## location = one of 4 pre-campaign locations an individual is from 


## Figure out which, if any,  demography elements will get added prior to the simulation 
## and what their affect will be on any of the subsequent models

## Figure out which variables will get assigned post simulation (mother recall/card indicator)
## These variables will depend on an individual's vaccination and exposure history

## For example: If vacc==1 then some proportion of those individuals have have a 1 for card indicator and mother recall
## We then feed this information into the model and have it recall these numbers

## James said we should simulate the full data before Keya manipulated and removed certain variables. 
## Does this mean we simulate the data set before momrec and docrec were combined into the overall recall variables?




#### Simulation Code

## Load necessary packages
devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)

##***************************1.1: Simulation Settings**************************** 
## Specify the number of time periods to simulate 
## Since children in Keya's data set are 0.7 to 14.99, let's simulate 15 years
## I decided to use monthly time steps since case study 1 model_pars parameters are 
## already structured for monthly time steps
times <- seq(1,180,by=1) 

## Set simulation settings needed for runserosim
simulation_settings <- list("t_start"=1,"t_end"=max(times))

## Specify the number of individuals in the simulation to match the number of individual's in Keya's data
N<-100
# N<-2570

##***************************1.2: Population Demography************************** 
## Generate the population demography tibble
## Specify options for each demography element and the distribution of each within the population
## We are interested in simulating a population where individuals have varying nutritional statuses and can reside in either of 2 locations
aux <- list("Group"=list("name"="group","options"=c("1", "2", "3", "4"), "distribution"=c(0.2638132,0.2245136,0.2587549,0.2529183)))


## Specify birth_times to match the age distribution of Keya's data 
## Create a vector with the number of individuals in each age class
# age.dist<-c(115,122,131,133,136,150,128,146,148,127,292,346,292,248,56)
# birth_distributiom<-function(sample_size){
#   birth_times<-NULL
#   for(age_dist in seq(sample_size)){
#     range_start<-(12*(age_dist-1))
#     range_end<-(12*(age_dist-1))+11
#     ss<-sample_size[age_dist]
#     birth_times_tmp <-sample(range_start:range_end,ss,replace = TRUE)
#     birth_times<-c(birth_times,birth_times_tmp)
#   }
#   birth_times
# }
# birth_times<-birth_distribution(age.dist)


## Let's assume that no individuals are removed from the population and set prob_removal to 0
demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=max(times), prob_removal=0, aux=aux)
# demography <- generate_pop_demography(N, times, birth_times, limit=0, removal_min=0, removal_max=max(times), prob_removal=0, aux=aux)
## Count the number of individuals in each group/location
demography %>% filter(times==1) %>% count(group) 


##********************************1.3: Biomarker Map*******************************
## Create biomarker map
## We are only interested in 3 exposure types (MCV1,MCV2 and natural infection) against one biomarker 
biomarker_map <- tibble(exposure_id=c(1,2,3),biomarker_id=c(1,1,1)) 


##***************************1.4: Force of Infection and Exposure Model***********
## Create an empty array to store the force of infection for all exposure types
foe_pars <- array(0, dim=c(n_distinct(demography$group),max(times),n_distinct(biomarker_map$exposure_id)))

## Specify the force of infection for exposure ID 1 which represents natural infection
foe_pars[,,1] <- 0.2

## Specify the force of vaccination for exposure ID 2 which represents MCV1 vaccination
foe_pars[,,2] <- 0.4

## Specify the force of vaccination for exposure ID 3 which represents MCV2 vaccination
foe_pars[,,3] <- 0.05

## I specified the same value for all time steps within foe_pars for simplicity but we can change to varying numbers to match real world settings. 

## Specify a simple exposure model which calculates the probability of exposure directly from the force of infection at that time step
## In this selected model, the probability of exposure is 1-exp(-FOI) where FOI is the force of infection at that time.
exposure_model<-exposure_model_simple_FOI


##********************************1.5: Immunity Model*****************************
## Specify immunity model within serosim function below 
immunity_model<-immunity_model_vacc_ifxn_titer_prot

## Specify which exposure IDs represent vaccination events 
vacc_exposures<-c(2,3)

## Specify the age at which an individual is eligible for MCV1 and MCV2 vaccination
vacc_age<-c(NA,9,12)

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
## DOUBLE CHECK WITH KEYA ##
max_vacc_events<-c(NA,1,1)

## Plot titer-mediated protection curve given parameters specified within model_pars for biomarker 1 which will be loaded in section 1.6
plot_titer_mediated_protection(0:7500, titer_prot_midpoint=5000, titer_prot_width=.001)
## These are the current parameters used within model_pars_cs1
## Maybe we should start of with simpler versions with high titer-mediated protection and no boosting events?


##****1.6: Antibody Model, Antibody Kinetics Parameters, and draw_parameters*****
## Specify the antibody model 
antibody_model<-antibody_model_biphasic

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model (Section 1.5), the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model (Section 1.7) are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
model_pars_path <- system.file("extdata", "model_pars_keya.csv", package = "serosim")
model_pars <- read.csv(file = model_pars_path, header = TRUE)
model_pars

## Specify the draw_parameters function to use 
draw_parameters<-draw_parameters_random_fx_titer_dep

## Plot titer dependent boosting effects given parameters specified within model_pars for biomarker 1 (measles)
plot_titer_dependent_boosting(start=0, end=1500, by=1, titer_ceiling_threshold=1000, titer_ceiling_gradient=0.0009)



##*****************1.7: Observation Model and observation_times*************
## Specify observation model to be used within runserosim 

## Specify the limits of detection for each biomarker for the continuous assays (lower detection limit is 8IU/I and the upper is 5000 IU/I)
bounds<-tibble(biomarker_id=c(1,1),name=c("lower_bound","upper_bound"),value=c(8,5000))


## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify observation_times (serological survey sampling design) to observe biomarker 1 (aka measles antibody titer) across all individuals at the end of the simulation
observation_times<-tibble(i=1:N,t=max(times), b=1)



##***************************1.9 Run Simulation*************************************
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
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
)


##***************************1.10 Generate Plots************************************
## Plot antibody states and exposure histories for 10 individuals 
plot_subset_individuals_history(res$antibody_states, res$exposure_histories_long, subset=5, demography)

## Plot exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)

## Plot exposure probabilities for all exposure types
plot_exposure_prob(res$exposure_probabilities_long)

## Plot antibody states for all individuals
plot_titers(res$antibody_states)

## Plot observed titers
plot_obs_titers_one_sample(res$observed_antibody_states)

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)


##*****************Generate Data Set with Keya's Format ***************************
##*## Post Simulation Variables
## mom_indicator:  binary indicator of if an individual had mother recall of any mcv dose among children with or without card
## card_indicator: binary indicator of if an individual had record (via card) of any mcv dose
## none_indicator = binary indicator of if an individual had neither card nor recall
## Combine runserosim outputs into a data set for Keya 
df<-res$observed_antibody_states %>% select(i,t,observed)
dem<-demography %>% filter(times==max(times)) %>% select(i,birth,group) %>% rename(location=group)
df<- right_join(df, dem, by="i")
df<- df %>% mutate(age=(t-birth)/12)
df <- df %>%  mutate(age_yr=floor(age)) %>% mutate(titer=ifelse(observed>0,log10(observed),0)) %>% select(i,t,age,age_yr,location, titer)


## Pull true vaccination and infection times from runserosim exposure history output
exp_hist <- res$exposure_histories_long %>%  filter(value==1)
inf_hist <- exp_hist %>% filter(x==1) %>% select(i,t) %>% rename(inf_time=t)
MCV1_hist <- exp_hist %>% filter(x==2) %>% select(i,t)  %>% rename(MCV1=t)
MCV2_hist <- exp_hist %>% filter(x==3) %>% select(i,t)  %>% rename(MCV2=t)

## Add vacc time to df 
df$MCV1<-NA
df$MCV2<-NA
df <- df %>% rows_update(MCV1_hist, by="i") %>% rows_update(MCV2_hist, by="i")

## Add days since vaccination
df<- df %>% mutate(daysbtwnmcvenr_1=180-MCV1, daysbtwnmcvenr_2=180-MCV2) 

## Add number of natrual infections to df (this will be hidden from Keya)
inf_count <- inf_hist %>% count(i) %>% rename(num_infs=n)
df$num_infs<-0
df<-df %>% rows_update(inf_count,by="i")
df_full<-df

## Individuals <=4 have some level vaccine information available 
df$vacc_obs <- ifelse(df$age_yr<=4,1,0)

## Add mom and card indicators to individuals with available vaccine information 
df$mom_indicator <- NA
df$card_indicator <- NA
for(indivs in df$i){
  if(df$vacc_obs[indivs]==0){
  df$daysbtwnmcvenr_1[indivs]<-NA
  df$daysbtwnmcvenr_2[indivs]<-NA
  }
  if(df$vacc_obs[indivs]==1){ ## If the individual is under 4 and we observed some level of vaccine information
    if(!is.na(df$MCV1[indivs])){  ## If the individual was actually vaccinated 
      df$mom_indicator[indivs] <- ifelse(runif(1)>0.7,1,0)  ## 30%  of individuals have accurate mother recall
      df$card_indicator[indivs] <- ifelse(runif(1)>0.8,1,0) ## 20%  of individuals have card info
    } else (df$card_indicator[indivs] <- ifelse(runif(1)>0.95,1,0)) ## If the individual was never actually vaccinated; false positive mother recall occurs 5% of the time
  }
}

## Add non indicator to individuals with available vaccine information 
df$none_indicator <- ifelse(df$mom_indicator==0 & df$card_indicator==0,1,0)

## Add boost_indicator
## 1 if individual received a vaccine within 1 year of the serosurvey
## 1 if individual<=2 with at least mcv1
## 1 if individual<=3 with mcv2
## 0 if none of the above characteristics are met 
df$boost_indicator<-ifelse(df$vacc_obs==0,0,
                           ifelse(df$daysbtwnmcvenr_1<=12,1,
                              ifelse(df$daysbtwnmcvenr_2<=12,1,
                                  ifelse(df$age<=2 & !is.na(df$daysbtwnmcvenr_1),1,
                                         ifelse(df$age<=3 & !is.na(df$daysbtwnmcvenr_2),1,0)))))
df$boost_indicator <-ifelse(is.na(df$boost_indicator),0,df$boost_indicator)
df_obs<-df %>% select(i,age,age_yr,location, titer, mom_indicator, card_indicator, none_indicator, boost_indicator)

df_full ## Full truth data set 
df_obs ## Observation data set (model input)
res ## All simulation outputs 


## Save outputs 
write.csv(df_full,"/Users/arthurmenezes/downloads//df_full.csv", row.names = FALSE)
write.csv(df_obs,"/Users/arthurmenezes/downloads//df_obs.csv", row.names = FALSE)
write.csv(res,"/Users/arthurmenezes/downloads//res.csv", row.names = FALSE)



## Graveyard 
# df$MCV1 <- NA
# df$MCV2 <- NA
# for(indivs in unique(vacc_hist$i)){
#   vacc_time_tmp<-vacc_hist$vacc_time[vacc_hist$i==indivs]
#   if(length(vacc_time_tmp)==1){
#     df$MCV1[indivs]<-vacc_time_tmp
#   }
#   if(length(vacc_time_tmp)==2){
#     df$MCV1[indivs]<-vacc_time_tmp[1]
#     df$MCV2[indivs]<-vacc_time_tmp[2]
#   }
#   if(length(vacc_time_tmp)==0){
#   }
# }
# 
# 
# ## Select a proportion of individuals within the serosurvey to be missing vaccine information 
# vacc_missing_prop<-0.7
# missing_indivs<-sample(1:N,floor(vacc_missing_prop*N),replace=FALSE)
# for(indivs in unique(missing_indivs)){
#   df$daysbtwnmcvenr_1[df$i==indivs]<-NA
#   df$daysbtwnmcvenr_2[df$i==indivs]<-NA
# }
# 
# 

