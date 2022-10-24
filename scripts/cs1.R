
devtools::load_all("~/Documents/GitHub/serosim")
devtools::document("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)

## Specify the number of time periods to simulate 
times <- seq(1,120,by=1) 

## Set simulation settings
simulation_settings <- list("t_start"=1,"t_end"=max(times))

N<-100

demography <- generate_pop_demography(N, times, limit=0, removal_min=0, removal_max=120, prob_removal=0)

antigen_map <- tibble(exposure_id=c(1,2),antigen_id=c(1,1)) 

## Create an empty array to store the force of infection for all exposure types
lambdas <- array(0, dim=c(1,max(times),n_distinct(antigen_map$exposure_id)))

## Specify the force of infection for exposure ID 1 which represents natural infection
lambdas[,,1] <- 0.2
lambdas[,,1] <- sample(seq(0,1,by=.01),max(times), replace=TRUE)

## Specify the force of vaccination for exposure ID 2 which represents vaccination
## Assume this is constant across all time steps
lambdas[,,2] <- 0.4
lambdas[,,2] <- sample(seq(0,1,by=.01),max(times), replace=TRUE)

## Specify a simple exposure model which calculates the probability of exposure directly from the force of infection at that time step
## The probability of exposure is 1-exp(-FOI) where FOI is the force of infection at that time
exposure_model<-exposure_model_simple_FOI

## Specify immunity model within serosim function below 
immunity_model<-immunity_model_vacc_ifxn_titer_prot

## Specify which exposure IDs represent vaccination events 
vacc_exposures<-2

## Specify the age at which an individual is eligible for vaccination (9 months old for measles)
vacc_age<-c(NA,9)

## Specify the maximum number of vaccines an individual can receive for each exposure types; note non vaccine exposures are listed as NAs
max_vacc_events<-c(NA,1)

## Specify the antibody model 
antibody_model<-antibody_model_biphasic

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model, the titer-ceiling parameters needed for draw_parameters and the observation error needed for the observation model are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
theta_path <- system.file("extdata", "theta_cs1.csv", package = "serosim")
theta <- read.csv(file = theta_path, header = TRUE)

## Specify the draw_parameters function to use 
draw_parameters<-draw_parameters_random_fx_titer_dep


## Limits of detection for continuous assays
boundary<-c(20,20000)

## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify observation_times to observe antigen 1 (aka Measles antibody titer) across all individuals at the midpoint and the end of the simulation (t=60 and t=120)
obs1 <- tibble(i=1:N,t=60, ag=1)
obs2 <- tibble(i=1:N,t=120, ag=1)
observation_times<-rbind(obs1,obs2)


res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  lambdas,
  antigen_map,
  theta,
  exposure_model,
  immunity_model,
  antibody_model,
  observation_model=observation_model_continuous_bounded_noise,
  draw_parameters=draw_parameters_random_fx_titer_dep,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
)


## Plot antibody states and exposure histories for 10 individuals 
plot_subset_individuals_history(res$antibody_states, res$exposure_histories_long, subset=6, demography)

## Plot exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)

## Plot exposure probabilities for all exposure types
plot_exposure_prob(res$exposure_probabilities_long)

## Plot antibody states for all individuals
plot_titers(res$antibody_states)

## Plot the first serosurvey at time 60 
obs60<-res$observed_antibody_states %>% filter(t==60)
plot_obs_titers_one_sample(obs60)

## Plot the second serosurvey at time 120 
obs120<-res$observed_antibody_states %>% filter(t==120)
plot_obs_titers_one_sample(obs120)

## Plot both serosurveys paired samples
plot_obs_titers_paired_sample(res$observed_antibody_states)

## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)


## Create plots for paper

library(cowplot)
plot_grid(plot_exposure_prob(res$exposure_probabilities_long), plot_exposure_histories(res$exposure_histories_long), nrow=1, ncol=2, align = "hv", scale=c(.98,.98))


plot_subset_individuals_history <- function(titers, exposure_histories, subset, demography){
  exposure_histories$e <- paste0("Exposure: ", exposure_histories$e)
  titers$ag <- paste0("Antigen: ", titers$ag)
  
  exposure_histories_subset<-exposure_histories %>% drop_na() %>% filter(value==1)
  removal_subset <- demography %>% filter(times==1)
  
  sample_indivs <- sample(1:N, size=subset)
  
  g<-  ggplot() +
    geom_vline(data=exposure_histories_subset %>% filter(i %in% sample_indivs), aes(xintercept=t, colour=e),linetype="dotted") +
    # geom_vline(data=removal_subset %>% filter(i %in% sample_indivs), aes(xintercept=removal, color="Removal Time"),linetype="solid") +
    geom_line(data=titers %>% filter(i %in% sample_indivs), aes(x=t,y=value,colour=ag)) +
    facet_wrap(~i) + theme_bw() +
    scale_color_hue("Key", guide=guide_legend(order=3)) +
    ggplot2::labs(title="Individual Antibody Kinetics",
                  x="Time",
                  y="Antibody Titer") + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
  return(g)
}

p3<-plot_subset_individuals_history(res$antibody_states, res$exposure_histories_long, subset=1, demography)
plot_grid(plot_titers(res$antibody_states), plot_obs_titers_one_sample(res$observed_antibody_states %>% filter(t==60)), p3, nrow=1, ncol=3, scale=c(.9,.9, .9))


plot_grid(plot_titers(res$antibody_states), plot_obs_titers_one_sample(res$observed_antibody_states %>% filter(t==60)), nrow=1, ncol=2, scale=c(.9,.9))



