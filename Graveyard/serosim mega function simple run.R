####LOAD NECESSARY PACKAGES#####
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(serosim)


#### ANTIBODY KINETICS PARAMATERS####
#Define antibody kinetics parameters; should be on the normal scale
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


#### CORE OF THE SIMULATION #####
completed_run<-run_full_simulation(500, 10, 12, 2, 120, 9, 60, kinetics_pars,
                                   overall_infection_probs=c(1, 0.75), titre_prot_midpoint=c(2000,4000), titre_prot_width=c(.01,.01), 0.7,
                                   pathogen, i, t)
attach(completed_run)



#### RESHAPE AND PLOT OUTPUTS####

#Go through the pre-calculated probabilities of infection matrices, reshape and plot
all_infection_probs <- reshape_all_infection_probs(all_infection_probs, 2)
plot_all_infection_probs(all_infection_probs)

#Go through the pre-calculated infection history matrices, reshape and plot
infection_histories <- reshape_infection_histories(infection_histories, 2)
plot_infection_histories(infection_histories)

#Go through the pre-calculated titre matrices, reshape and plot
titres<- reshape_titres(titres, 2)
plot_titres(titres)
