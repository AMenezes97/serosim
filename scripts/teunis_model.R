setwd("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
devtools::document()
devtools::load_all()

## Let's simulate 1 pathogen (pertussis) with the Teunis-style antibody model
## We'll generate the FOI from an SIR model

## 1. Simulation settings
times <- seq(1,365*2,by=1) ## 2 years in daily resolution
simulation_settings <- list("t_start"=min(times),"t_end"=max(times))

## 2. Population demography
N <- 100
birth_times <- simulate_birth_times(N, times, age_min=0)
removal_times <- simulate_removal_times(N, times=times,birth_times=birth_times, removal_min=365,removal_max=max(times),prob_removal=1)

aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "proportion"=c(0.5,0.5)),
            "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "proportion"=c(0.25,0.25,0.25,0.25)),
            "SES"=list("name"="ses","options"=c("low","high"),"proportion"=c(0.8,0.2)))
#demography <- generate_pop_demography(N=N,times=times, age_min=180,prob_removal=1,aux=aux)
demography <- generate_pop_demography(N=N,times=times, age_min=180,prob_removal=0)

biomarker_map <- tibble(exposure_id=c("infection"),biomarker_id=c("IgG"))
biomarker_map <- reformat_biomarker_map(biomarker_map)

## 3. FOE models
foe_pars <- data.frame(x=1,g=1,name=c("beta","gamma","I0","R0","t0"),values=c(0.2,1/7,1/10000,0,50))
plot_exposure_model(exposure_model_sir,seq(1,365,by=1),n_groups=1,n_exposures=1,foe_pars=foe_pars)

## 4. Biomarker map
biomarker_map <- tibble(exposure_id=c(1),biomarker_id=c(1))

## 5. Immunity model
immunity_model_all_successful

## 6. Antibody model
model_pars <- read_csv("inst/extdata/model_pars_typhoid.csv")
#plot_antibody_model(antibody_model_typhoid,N = 25,times=seq(1,365,by=1),model_pars=model_pars,
#                    biomarker_map=biomarker_map,draw_parameters_fn=draw_parameters_random_fx)

bounds<-data.frame(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(10,500))

## Run serosim
set.seed(1)
res<- runserosim(
    simulation_settings,
    demography, 
    observation_times=NULL,
    foe_pars, 
    biomarker_map,
    as.data.frame(model_pars),
    exposure_model=exposure_model_sir, 
    immunity_model=immunity_model_all_successful, 
    antibody_model=antibody_model_typhoid, 
    observation_model=observation_model_continuous_bounded_noise,
    draw_parameters=draw_parameters_random_fx, 
    
    ## Pre-specified parameters/events
    exposure_histories_fixed=NULL,
    bounds=bounds,
    VERBOSE=TRUE
)
plot_exposure_histories(res$exposure_histories_long)
plot_titers(res$observed_antibody_states)
 
res$observed_antibody_states %>% filter(t == 730) %>% ggplot() + geom_jitter(aes(x=as.factor(t),y=observed),width=0.25)
