####LOAD NECESSARY PACKAGES#####
# library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

####LOAD PACKAGE#####
library(serosim) #serosim
devtools::load_all("../../serosim")


#### POPULATION PARAMETERS####
N <- 10 # number of individuals in the study population
years <- 1 # number of years to run the simulation for
t_periods_per_year <- 12 # number of time periods per year (12 months)
times<- seq(1,years*t_periods_per_year, by=1) # creates a vector of all months in 10 years
N_pathogens<-2
obs_time <-120 # time step when the serosurvey is conducted 
vacc_min <- 9 # time step when an individual becomes eligible for vaccination 
vacc_max <- 60 #time step at which an individual becomes too old to get vaccinated 

birth_times <-simulate_birth_times(N, times, limit=vacc_min)
# ?simulate_birth_times
# N_alive <- find_N_alive(N, times)  #This function doesn't work when individual's are removed from the population 
# ?find_N_alive
plot_age_distribution(birth_times)
# ?plot_age_distribution



removal_min <- 0 # the minimum age at which an individual can be removed from the population 
removal_max <- 100 # the maximum age at which an individual can be removed from the population 
prob_removal <- 0.9 # the probability that an individual will be removed from the population

removal_times<-simulate_removal_times(N, times, birth_times, removal_min=0, removal_max=100, prob_removal=0.9)
# ?simulate_removal_times
plot_age_distribution(removal_times)


#pre load the demography categories, values and distributions 
aux <- list("SES"=list("name"="SES","options"=c("low","medium","high"), "distribution"=c(0.2,0.2,0.6)),
            "NS"=list("name"="NS","options"=c("low","high"),"distribution"=c(0.5,0.5)),
            "Sex"=list("name"="Sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),
            "Location"=list("name"="Location","options"=c("North", "South", "East", "West"), "distribution"=c(0.25,0.25,0.25,0.25)) )

demography<- generate_pop_demography(N, times, birth_times, removal_times, aux=aux)



#### ANTIBODY KINETICS PARAMATERS####
#Define antibody kinetics parameters 
#Antibody kinetics parameters should be on the normal scale
# ?load_kinetics_parameters
kinetics_pars<-load_kinetics_parameters(boost_infection_long_mean_1=9000, boost_infection_long_var_1=500,
                                         wane_infection_long_mean_1=0, wane_infection_long_var_1=0,
                                         boost_infection_short_mean_1=3000, boost_infection_short_var_1=500,
                                         wane_infection_short_mean_1=1/120, wane_infection_short_var_1=0,
                                         boost_vacc_long_mean_1=3000, boost_vacc_long_var_1=500,
                                         wane_vacc_long_mean_1=0, wane_vacc_long_var_1=0,
                                         boost_vacc_short_mean_1=1000,boost_vacc_short_var_1=1000,
                                         wane_vacc_short_mean_1=1/120, wane_vacc_short_var_1=0,
                                         titre_ceiling_gradient_1=0.5/200, titre_ceiling_threshold_1=200,
                                         boost_infection_long_mean_2=9000,boost_infection_long_var_2=500,
                                         wane_infection_long_mean_2=0, wane_infection_long_var_2=0,
                                         boost_infection_short_mean_2=3000, boost_infection_short_var_2=500,
                                         wane_infection_short_mean_2=1/120, wane_infection_short_var_2=0,
                                         boost_vacc_long_mean_2=3000, boost_vacc_long_var_2=500,
                                         wane_vacc_long_mean_2=0,wane_vacc_long_var_2=0,
                                         boost_vacc_short_mean_2=1000, boost_vacc_short_var_2=500,
                                         wane_vacc_short_mean_2=1/120, wane_vacc_short_var_2=0,
                                         titre_ceiling_gradient_2=0.5/200, titre_ceiling_threshold_2=200)


# ?plot_titre_dependent_boosting
#Relationship between titre and boosting for Pathogen 1
plot_titre_dependent_boosting(0,1000, by=10, kinetics_pars = kinetics_pars, pathogen=1)
#Relationship between titre and boosting for Pathogen 2
plot_titre_dependent_boosting(0,1000, by=10, kinetics_pars = kinetics_pars, pathogen=2)




#### FORCE OF INFECTION####
#For two pathogens, assume that cumulative probability of infection is 100% and 75% over the
#period of the simulation.
overall_infection_probs <- c(1, 0.75) 


#Option 1: Random force of infection for all pathogens 
FOIs <- draw_samples(times, N_pathogens, kernel_fn = se_kernel, sigma=2,length = 24)
# ?generate_prob_infections
prob_infections<-generate_prob_infections(FOIs, N_pathogens, overall_infection_probs)
# ?plot_FOIs
plot_FOIs(prob_infections)
  
#Option 2: Constant circulation of pathogen 1 and no circulation of pathogen 2 
# prob_infections<-generate_prob_infections(pathogens=N_pathogens, overall_infection_probs=overall_infection_probs, type=2, times=times, set_prob1=0.07, set_prob2 = 0)
# plot_FOIs(prob_infections)

#Option 3: Constant circulation of both pathogens (note that the pathogens have different overall probabilities of infection so probabilities of infections will ultimately be different)
# prob_infections<-generate_prob_infections(pathogens=N_pathogens, overall_infection_probs=overall_infection_probs, type=2, times=times,  set_prob1=0.09, set_prob2 = 0.3)
# plot_FOIs(prob_infections)

#Option 4: Input your own FOI from SIR model 


#### TITRE-MEDIATED PROTECTION####
#Parameters of the titre protection curve.
titre_prot_midpoint<- c(2000,4000) #You have 50% protection at a titre of 2000 for pathogen 1 and 50% protection at a titre of 2000 for pathogen 2
titre_prot_width <- c(.01,.01) 
# ?plot_titre_mediated_protection
plot_titre_mediated_protection(0:10000, N_pathogens, titre_prot_midpoint, titre_prot_width)

#### VACCINATION HISTORIES####

#The overall probability of vaccination
prob_vaccination <- 0.7

#Simulate vaccine histories
# ?simulate_vaccine_histories
vaccine_histories<-simulate_vaccine_histories(N, times, birth_times, vacc_min, vacc_max, prob_vaccination)

#Plot vaccine histories
# ?plot_vaccine_histories
plot_vaccine_histories(vaccine_histories)


#### CORE OF THE SIMULATION #####
##********************************************

completed_matrices<- run_simulation1(birth_times, N, times, N_pathogens, kinetics_pars, vaccine_histories, pathogen, i, t)
attach(completed_matrices)


#### SAVE FULL TRUTH OUTPUTS#### 
all_infection_probs_FT <- all_infection_probs
infection_histories_FT <- infection_histories
titres_FT <- titres


#### RESHAPE AND PLOT OUTPUTS####

#Go through the pre-calculated probabilities of infection matrices, reshape and plot 
all_infection_probs <- reshape_all_infection_probs(all_infection_probs, N_pathogens)
plot_all_infection_probs(all_infection_probs)

#Go through the pre-calculated infection history matrices, reshape and plot 
infection_histories <- reshape_infection_histories(infection_histories, N_pathogens)
plot_infection_histories(infection_histories)

#Go through the pre-calculated titre matrices, reshape and plot 
titres<- reshape_titres(titres, N_pathogens)
plot_titres(titres)

#Reshape age matrix for use later on 
age<-reshape_age(age, N_pathogens)

#Reshape vaccine history matrix for plotting 
vaccine_histories_reshaped <- reshape_vaccine_histories(vaccine_histories)

##### OBSERVATION PROCESS####
#Add a column to the titres data set to simulate the observation process and any noise attributed to assay variability
obs_sd <- 0.25
titres<-observation_process(titres, obs_sd)


#### GENERATE PLOTS####

#Plot the infection and vaccination history for a subset of individuals
# ?plot_subset_individuals_history
plot_subset_individuals_history(infection_histories, vaccine_histories_reshaped, titres, 10, N_pathogens)

#Plot serosurvey results 
# ? plot_serosurvey_observation
plot_serosurvey_observation(titres, obs_time, titre="true") 
plot_serosurvey_observation(titres, obs_time, titre="observed")

#Generate a serosurvey data set with tires at sampling time and pull relevant information from the simulation outputs 
# ? full_truth
titresobs<- full_truth(titres, obs_time, infection_histories, vaccine_histories_reshaped, age, N_pathogens)


#Plot serosurvey results and label by true events
# ? plot_serosurvey_observation
plot_serosurvey_observation(titresobs, obs_time, titre="true", Events="true") 
plot_serosurvey_observation(titresobs, obs_time, titre="observed", Events="true")


# ?plot_serosurvey_distribution
plot_serosurvey_distribution(titresobs, titre="true", seronegative="keep", bin= 200, color="pathogen" )
plot_serosurvey_distribution(titresobs, titre="true", seronegative="remove", bin= 200, color="pathogen" )
plot_serosurvey_distribution(titresobs, titre="true", seronegative="keep", bin= 200, color="events" )
plot_serosurvey_distribution(titresobs, titre="true", seronegative="remove", bin= 200, color="events")

