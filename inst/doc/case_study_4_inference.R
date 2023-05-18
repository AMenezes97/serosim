## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7,fig.height=5
)

## ----message=FALSE, warning=FALSE, echo=FALSE,eval=TRUE-----------------------
## Install and load serosim 
## devtools::install_github("seroanalytics/serosim")
library(serosim)

## Load required packages 
library(plyr) ## Needed for a few internal serosolver functions -- careful, as some conflicts with dplyr
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(coda)
library(foreach)
library(doParallel)
library(parallel)
library(viridis)
library(bayesplot)
library(ggpubr)

## Custom inference packages -- you will need a C++ tool chain i.e. Rtools
## Install serofoi
## Make sure rstan is installed
library(rstan)
## devtools::install_github("epiverse-trace/serofoi")
library(remotes)
remotes::install_github("epiverse-trace/serofoi")
library(serofoi)

## Install serosolver
remotes::install_github("seroanalytics/serosolver")
library(serosolver)

set.seed(123)

## -----------------------------------------------------------------------------
## Specify the number of time periods to simulate 
times <- seq(1,50,by=1) 
simulation_settings <- list("t_start"=1,"t_end"=max(times))

## Generate the population demography tibble
## Specify the number of individuals in the simulation; N=1000
N <- 1000

## Assign individuals into either urban or rural location
aux_demography <- list("Group"=list("name"="location","options"=c("Urban","Rural"),"distribution"=c(0.7,0.3)))
demography <- generate_pop_demography(N=N, times=times,aux=aux_demography,age_min=0,prob_removal=0)
demography$group <- as.numeric(factor(demography$location, levels=c("Urban","Rural"))) ## Convert the location ID to the group ID for the main simulation.

## Vaccination will be determined by age and doses, not dependent on biomarker quantity
max_events <- c(10,2) ## Maximum of 10 exposure ID 1 events (infection), and 2 exposure ID 2 events (vaccination).
vacc_exposures <- 2 ## Vaccination is the 2nd exposure ID
vacc_age <- c(NA,2) ## Individuals can be infected at any age (vacc_age[1] = NA), but only vaccinated once they reach 2 years old (vacc_age[2]=2)

# Create simple biomarker map
biomarker_map <- tibble(exposure_id=c("ifxn","vacc"),biomarker_id=c("IgG","IgG"))
## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map)


## -----------------------------------------------------------------------------
## Create an empty array to store the force of exposure for all exposure types
# Dimension 1: location
# Dimension 2: time
# Dimension 3: exposure ID
foe_pars <- array(0, dim=c(2,max(times),n_distinct(biomarker_map$exposure_id)))

## Specify the force of exposure for Location 1 (urban), Exposure ID 1 which represents natural infection
foe_pars[1,,1] <- c(rep(0.05,20),rep(0.025,15),rep(0.01,15))
#foe_pars[1,,1] <-0.03

## Specify the force of exposure for Location 2 (rural), Exposure ID 1 which represents natural infection
foe_pars[2,,1] <- rev(c(rep(0.05,20),rep(0.025,15),rep(0.01,15))*2)
#foe_pars[2,,1] <- 0.06

## Specify a simple exposure model which calculates the probability of exposure 
## directly from the force of exposure at that time step. In this selected model, 
## the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.
exposure_model<-exposure_model_simple_FOE

## Examine the probability of exposure over time for the specified exposure model and foe_pars array
plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times, 
                    n_groups = 2, n_exposures = 2, foe_pars=foe_pars)

## -----------------------------------------------------------------------------
## Assume that immunity from infection is dependent on latent biomarker level at time of exposure
p_immunity <- plot_biomarker_mediated_protection(seq(0,8,by=0.1), 2, 2)

## Tibble of model parameter controls related to infection immunity
model_pars_immunity <- tibble("exposure_id"="ifxn","biomarker_id"="IgG",
                              "name"=c("biomarker_prot_midpoint","biomarker_prot_width"),
                              "mean"=c(2,2),"sd"=NA,"distribution"="")
## See example in help file
immunity_model <- immunity_model_vacc_ifxn_biomarker_prot

## Specify the antibody model 
antibody_model<-antibody_model_monophasic

## Bring in the antibody parameters needed for the antibody model
model_pars_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars <- reformat_biomarker_map(bind_rows(model_pars_original,model_pars_immunity))

model_pars[model_pars$name == "wane" & model_pars$exposure_id == 1,"mean"] <- 0.008
model_pars[model_pars$name == "wane" & model_pars$exposure_id == 2,"mean"] <- 0.005
model_pars[model_pars$name == "wane","sd"] <- 0.025
## Specify the draw_parameters function
draw_parameters<-draw_parameters_random_fx

## Plot example biomarker trajectories given the specified antibody kinetics model, 
## model parameters and draw parameters function 
p_antibody <- plot_antibody_model(antibody_model_monophasic, N=100, model_pars=model_pars %>% drop_na(),
                    draw_parameters_fn = draw_parameters_random_fx, 
                    biomarker_map=biomarker_map)
p_immunity + ggtitle("Relationship between pre-exposure biomarker level\n and probability of infection")
p_antibody + ggtitle("100 randomly drawn post-exposure biomarker trajectories")


## -----------------------------------------------------------------------------
## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise

## Specify assay sensitivity and specificity needed for the observation model
model_pars_original[model_pars_original$name =="obs_sd","sd"] <- 0.5
sensitivity<-0.85
specificity<-0.95

bounds <- dplyr::tibble(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(1,10))

## Specify observation_times (serological survey sampling design) to observe
## biomarker 1 at two time points around t 80 and 110
observation_times <- tibble(i=rep(1:max(demography$i),each=2),t=rep(c(48,50),N))

## -----------------------------------------------------------------------------
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
  bounds=bounds,
  max_events=max_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  sensitivity=sensitivity,
  specificity=specificity
)

## ----message=FALSE, warning=FALSE, results=FALSE------------------------------
library(serofoi)
## Pull our the observed biomarker states and combine with the demography data to get each individual's age
serodat <- res$observed_biomarker_states %>% select(i, t, observed) %>% filter(t == 50)
demog <- res$demography %>% select(i, birth, location) %>% distinct()
serodat <- left_join(serodat, demog, by="i") 

## Create a custom function for summarizing the simulated serosurvey data into seroprevalence by age,
## then fit the serocatalytic model to the Rural and Urban location data.
fit_serofoi <- function(seropos_cutoff = 1.5, serodat, true_foi){
  summary <- serodat %>% dplyr::mutate(age = t - birth) %>% 
    dplyr::mutate(seropos=observed >= seropos_cutoff) %>% 
    group_by(age, location, t) %>% 
    dplyr::summarise(total=n(),counts=sum(seropos),.groups = "drop") %>% ungroup() %>%
    dplyr::mutate(age_min = age, age_max=age) %>%
    dplyr::mutate(test="sim",antibody="IgG",country="sim") %>%
    dplyr::rename(tsur=t, survey=location)
  
  serodat_rural <- prepare_serodata(summary %>% filter(survey=="Rural"))
  model_rural <- run_seromodel(serodat_rural,foi_model = "tv_normal",n_iters=3000)
  
  p_rural <- plot_seromodel(model_rural, size_text = 8,foi_sim=true_foi[2,,1])
  
  serodat_urban <- prepare_serodata(summary %>% filter(survey=="Urban"))
  model_urban <- run_seromodel(serodat_urban,foi_model = "tv_normal",n_iters=3000)
  p_urban <- plot_seromodel(model_urban, size_text = 8,foi_sim=true_foi[1,,1])
  p_rural | p_urban
}

## Compare results for different seropositivity thresholds
p1 <- fit_serofoi(1.25,serodat, foe_pars)
p2 <- fit_serofoi(2,serodat, foe_pars)

## Seropositivity threshold of 1.25
p1

## Seropositivity threshold of 2
p2

## ---- eval=FALSE--------------------------------------------------------------
#  ## Make the antibody waning rates one order of magnitude slower with little variation
#  model_pars[model_pars$name == "wane" & model_pars$exposure_id == 1,"mean"] <- 0.0008
#  model_pars[model_pars$name == "wane" & model_pars$exposure_id == 2,"mean"] <- 0.0005
#  model_pars[model_pars$name == "wane","sd"] <- 0.00025

## -----------------------------------------------------------------------------
## Specify the force of exposure for exposure ID 2 which represents vaccination
foe_pars[1,,2] <- 0.25
foe_pars[2,,2] <- 0.01

plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times, 
                    n_groups = 2, n_exposures = 2, foe_pars=foe_pars)


## ----message=FALSE, warning=FALSE, results=FALSE------------------------------
## Run the core simulation and save outputs in "res"
res_vaccination <- runserosim(simulation_settings,demography,observation_times,
  foe_pars,biomarker_map,model_pars,
  exposure_model,immunity_model,antibody_model,observation_model,
  draw_parameters,bounds=bounds,max_events=max_events,
  vacc_exposures=vacc_exposures,vacc_age=vacc_age,
  sensitivity=sensitivity,specificity=specificity)

serodat_vaccine <- res_vaccination$observed_biomarker_states %>% select(i, t, observed) %>% filter(t == 50)
demog <- res_vaccination$demography %>% select(i, birth, location) %>% distinct()
serodat_vaccine <- left_join(serodat_vaccine, demog, by="i") 
p_vaccine_foi <- fit_serofoi(2,serodat_vaccine, foe_pars)
p_vaccine_foi

## -----------------------------------------------------------------------------
library(serosolver)

run_name <- "serosim_recovery" ## Name to give to created files
main_wd <- getwd()

## Create a directory to store the MCMC chains
chain_wd <- paste0(main_wd,"/chains/",run_name)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

#options(mc.cores=5)
n_chains <- 5 ## Number of MCMC chains to run

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", chain_wd))

## Run multiple chains in parallel
## R CMD check only allows a maximum of two cores 
cl <- makeCluster(2)
registerDoParallel(cl)

## ---- echo=FALSE--------------------------------------------------------------
## Create some custom functions to use for this vignette
## Converts serosim exposure history into the format expected by serosolver
## Uses strain_isolation_times to bucket exposures into time blocks
convert_inf_hist_to_serosolver <- function(true_inf_hist, strain_isolation_times){
  true_inf_hist$t_group <- cut(true_inf_hist$t, breaks=c(strain_isolation_times,max(strain_isolation_times)+1),include.lowest=TRUE,right=FALSE)
  true_inf_hist <- true_inf_hist %>% group_by(t_group) %>% dplyr::mutate(t_floor = min(t))
  true_inf_hist <- true_inf_hist %>% mutate(value = ifelse(is.na(value), 0, value)) %>%
    dplyr::select(-t) %>%
    dplyr::rename(t=t_floor) %>%
    group_by(i, t) %>% dplyr::summarize(x=sum(value)) %>% mutate(x = ifelse(x >=1, 1, x))
  true_inf_hist <- as.matrix(acast(true_inf_hist, formula=i ~ t))
  true_inf_hist
}
## Converts serosim sero data into the format expected by serosolver
## time_blocks is used if we want to run serosolver at a lower time resolution than the serosim model
convert_serodata_to_serosolver <- function(sero_data, time_blocks=1){
    add_val <- 0
  if(time_blocks > 1) add_val <- 1
  colnames(sero_data) <- c("individual","samples","virus","true_titre","titre","DOB","removal","sex","location")
  sero_data$run <- 1 ## Variable used if we have multiple observations per time
  sero_data$group <- 1 ## Variable used if we want to stratify the population into distinct groups
  sero_data <- sero_data %>% select(individual,samples,virus,titre,DOB,run,group)
  sero_data$samples <- floor(sero_data$samples/time_blocks) + add_val
  sero_data$DOB <- floor(sero_data$DOB/time_blocks)  + add_val
  as.data.frame(sero_data)
}
## Expand the sero data to allow model to be solved at all times
expand_sero_data <- function(sero_data, strain_isolation_times){
  sero_data_tmp <- expand_grid(individual=unique(sero_data$individual),
                               samples=strain_isolation_times,virus=1,group=1,run=1)
  sero_data_tmp <- sero_data_tmp %>% 
    left_join(sero_data %>% select(individual, DOB) %>% distinct()) %>% 
    left_join(sero_data %>% select(individual, samples,titre) %>% distinct()) %>% 
    filter(samples >= DOB)
  sero_data_tmp
}

## -----------------------------------------------------------------------------
## Some data cleaning to convert the serosim output into the format expected by serosolver
sero_data <- convert_serodata_to_serosolver(res$observed_biomarker_states %>% left_join(demography %>% select(-times) %>% distinct(),by="i"))
## Note each individual's group/location
sero_data <- sero_data %>% select(-group) %>% left_join(res$demography %>% select(i, group) %>% distinct() %>% dplyr::rename(individual=i), by="individual")

## -----------------------------------------------------------------------------
## Set up parameter table
par_tab <- read.csv("par_tab_base.csv",stringsAsFactors=FALSE)

## -----------------------------------------------------------------------------
## These are the alpha and beta priors on the Beta distribution prior on the per-time attack rate. 
## Set to 1/1 for uniform, can e.g., set to 1/100 to represent a strong prior on low per-time attack rate.
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(1,1)

## -----------------------------------------------------------------------------
## Setup some inputs for serosolver -- this just tells the function the vector of possible
## exposure times. The antigenic map is uninformative here, but is used in other examples
## to capture cross-reactivity
#exposure_times <- seq(1,max(sero_data$samples),by=1) # If we wanted to estimated exposure states for each of 50 years
exposure_times <- c(seq(1,45,by=5), 46:50)
true_inf_hist <- res$exposure_histories_long
true_inf_hist$t_group <- cut(true_inf_hist$t, breaks=exposure_times) # Note which time block each exposure time is in

## THis is a data structure used by serosolver for antigenically variable pathogens -- it can be largely ignored here
antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=exposure_times)

## Total number of exposures
sum(true_inf_hist$value,na.rm=TRUE)

true_inf_hist <- convert_inf_hist_to_serosolver(true_inf_hist,exposure_times) # Convert the true exposure history to the format expected by serosolver

## -----------------------------------------------------------------------------
prior_func <- function(par_tab){
  f <- function(pars){
    pr <- dlnorm(pars["wane"],log(0.008),0.25,log=TRUE)
    pr
  }
}


## ---- echo=FALSE--------------------------------------------------------------
## MCMC settings, not super important but can be tweaked.
## It's fairly involved tweaking these to optimize the MCMC chains, but happy to chat through
## Main one to change is iterations
mcmc_pars <- c("save_block"=100,"thin"=20,"thin_hist"=200,
               "iterations"=500000,
               "adaptive_period"=500000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.5,
               "year_swap_propn"=1,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=2, "hist_opt"=1,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

## -----------------------------------------------------------------------------
## Set up posterior function for later. Data_type=2 is for continuous data
f <- create_posterior_func(par_tab,sero_data,version=2,data_type=2,antigenic_map=antigenic_map)
## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
output <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","serosolver")) %dopar% {
    index <- 1
    ## Try generating random starting values until we find a set that returns a finite likelihood
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
        start_tab <- generate_start_tab(par_tab)
        start_inf <- setup_infection_histories_total(sero_data,exposure_times,1,100)
        
        inf_hist_correct <- sum(check_inf_hist(sero_data, exposure_times, start_inf))
        y <- f(start_tab$values, start_inf)
        lik <- sum(y[[1]])
        index <- index + 1
    }
    ## Run serosolver!
    output <- serosolver::run_MCMC(start_tab, sero_data,
                                strain_isolation_times = exposure_times,
                                start_inf_hist=start_inf,
                                filename=x,
                                antigenic_map=antigenic_map,
                                CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                CREATE_PRIOR_FUNC = prior_func,
                                version=2,
                                mcmc_pars=mcmc_pars,
                                data_type=2)
}


## -----------------------------------------------------------------------------
## Read in the MCMC chains and create a trace plot
chains <- load_mcmc_chains(chain_wd,par_tab=par_tab,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
## Save traceplots from coda package
chains$theta_chain %>% as_tibble() %>% 
  dplyr::select(sampno, chain_no, mu_short, wane, error,total_infections) %>%
  pivot_longer(-c(sampno,chain_no)) %>%
  mutate(chain_no=as.factor(chain_no)) %>%
  ggplot() + geom_line(aes(x=sampno,y=value,col=chain_no)) +
  facet_wrap(~name,scales="free_y",ncol=2)

## Check Rhat statistic for antibody kinetics parameters
list_chains <- chains$theta_list_chains
list_chains <- lapply(list_chains, function(x) x[,(colnames(x) %in% c("mu_short","wane","error","total_infections"))])
gelman.diag(as.mcmc.list(list_chains))
effectiveSize(as.mcmc.list(list_chains))

## ---- warning=FALSE, message=FALSE,results=FALSE------------------------------
## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 
## Total number of infections

## ---- warning=FALSE, message=FALSE--------------------------------------------
print(paste0("True number of exposures: ", sum(true_inf_hist)))
print(paste0("Estimated median and 95% CrI on number of exposures: ",
             median(chains$theta_chain$total_infections), " (",
             quantile(chains$theta_chain$total_infections, 0.025), "-",
             quantile(chains$theta_chain$total_infections, 0.975), ")"))

print(paste0("True boosting mean: ", model_pars[model_pars$exposure_id == 1 & model_pars$name == "boost", "mean"]))
print(paste0("Estimated boosting (posterior mean and 95% CrI): ",  signif(mean(chains$theta_chain$mu_short),3), " (",
             signif(quantile(chains$theta_chain$mu_short, 0.025),3), "-",
             signif(quantile(chains$theta_chain$mu_short, 0.975),3), ")"))

print(paste0("True waning rate mean: ", model_pars[model_pars$exposure_id == 1 & model_pars$name == "wane", "mean"]))
print(paste0("Estimated waning rate (posterior mean and 95% CrI): ",  signif(mean(chains$theta_chain$wane),3), " (",
             signif(quantile(chains$theta_chain$wane, 0.025),3), "-",
             signif(quantile(chains$theta_chain$wane, 0.975),3), ")"))


## -----------------------------------------------------------------------------
## Plot individual infection history estimates
p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:25,
                                             real_inf_hist = as.matrix(true_inf_hist),
                                             strain_isolation_times = exposure_times,nsamp=100,
                                             number_col = 5)
p_cumu_infs[[1]]

## -----------------------------------------------------------------------------
p_cumu_infs[[2]]

## -----------------------------------------------------------------------------
## Plot attack rates
## Get number alive in each time point and true number of exposures per time point
inf_chain <- inf_chain %>% left_join(res$demography %>% select(i, group) %>% distinct())

## Get number of individuals alive to be exposed in each time period by group
n_alive <- get_n_alive_group(sero_data, exposure_times)
colnames(true_inf_hist) <- exposure_times
## Get true number of exposures from serosim output by time period
n_inf <- true_inf_hist %>% 
  as_tibble() %>% 
  dplyr::mutate(i = 1:n()) %>% 
  pivot_longer(-i) %>% 
  dplyr::mutate(name = as.numeric(name)) %>%
  left_join(res$demography %>% select(i, group) %>% distinct()) %>% group_by(name, group) %>% 
  dplyr::summarize(ar=sum(value))

n_alive <- n_alive %>% 
  as_tibble() %>% 
  dplyr::mutate(group=1:n()) %>% 
  pivot_longer(-group) %>% 
  dplyr::mutate(name=as.numeric(name))

## Create tibble of true per-time attack rates
true_ar <- n_inf %>% left_join(n_alive) %>% mutate(AR = ar/value) %>% dplyr::rename(j=name)

## Group 2 is the rural location, group 1 is the urban location
p_ar <- plot_attack_rates(inf_chain, sero_data, exposure_times,true_ar=true_ar,by_group=TRUE,pad_chain=FALSE)
p_ar

## ----fig.height=7,fig.width=8-------------------------------------------------
## Plot model fits to titre data. We expand the sero_data object so that the function plots
## for all possible times, not just those with observation times.
sero_data_tmp <- expand_sero_data(sero_data, exposure_times)

titre_pred_p <- plot_infection_histories(chain = chain[chain$chain_no == 1,], 
                         infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                         titre_dat = sero_data_tmp, 
                         individuals = 1:25,
                         strain_isolation_times = exposure_times,
                         nsamp = 100, # Needs to be smaller than length of sampled chain 
                         par_tab = par_tab,p_ncol=5, data_type=2)
titre_pred_p

