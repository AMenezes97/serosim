README
================

# Motivation

The serosim package is designed to simulate serological survey data
arising from user-specified vaccine or infection-generated and antibody
kinetics processes. serosim allows users to specify and adjust model
inputs responsible for generating the observed titer values like
time-varying patterns of infection and vaccination, population
demography, immunity and antibody kinetics, and serological survey
sampling design in order to best represent the population and disease
system(s) of interest.

Here, we will use the serosim package to generate a simple
cross-sectional serosurvey at the end of a 10 year simulation period for
100 individuals who have either been vaccinated, infected, or both. We
will set up each of the required arguments and models for the runserosim
function in the order outlined in the methods section of the paper. We
then run the simulation and examine its outputs

# Installation

Load necessary packages:

``` r
devtools::load_all("~/Documents/GitHub/serosim")
library(tidyverse)
library(data.table)
library(ggplot2)
```

# 1.1 Simulation Settings

We will simulate monthly time steps across a 10 year period. Therefore,
we will have 120 time steps. Note that these are arbitrary time steps
which will need to be scaled to the right time resolution to match any
time-based parameters used later on.

``` r
## Specify the number of time steps in the simulation
times <- seq(1,120,by=1) 

## Set simulation settings argument needed for runserosim 
simulation_settings <- list("t_start"=1,"t_end"=max(times))
```

# 1.2 Population Demography

For this case study, we are not interested in tracking any demography
information other than an individual’s birth time. We will use the
generate_pop_demography function to create the demography tibble needed
for runserosim.

Note: The runserosim function only requires a demography tibble with two
columns (individuals and times). In this case it will assume all
individuals are born at the start of the simulation period.

``` r
## Generate the population demography tibble
## Specify the number of individuals in the simulation; N=100
## No individuals are removed from the population; prob_removal=0
demography <- generate_pop_demography(N=100, times=times, removal_min=0, removal_max=120, prob_removal=0)

## Examine the generated demography tibble
summary(demography)
```

    ##        i              times            birth           removal     
    ##  Min.   :  1.00   Min.   :  1.00   Min.   :  1.00   Min.   : NA    
    ##  1st Qu.: 25.75   1st Qu.: 30.75   1st Qu.: 29.00   1st Qu.: NA    
    ##  Median : 50.50   Median : 60.50   Median : 56.00   Median : NA    
    ##  Mean   : 50.50   Mean   : 60.50   Mean   : 57.53   Mean   :NaN    
    ##  3rd Qu.: 75.25   3rd Qu.: 90.25   3rd Qu.: 89.00   3rd Qu.: NA    
    ##  Max.   :100.00   Max.   :120.00   Max.   :119.00   Max.   : NA    
    ##                                                     NA's   :12000

# 1.3 Antigen Map

Set up the exposure IDs and antigen IDs for the simulation which will
determine which infection or vaccination events are occurring. Here, we
will simulate one circulating pathogen (exposure_ID=1) and one vaccine
(exposure_ID=2)both of which will boost titers to the same
antigen(antigen_ID=1). This antigen map can be used for any simulations
of vaccine preventable diseases like measles vaccination and measles
natural infection.

``` r
## Create antigen map
antigen_map <- tibble(exposure_id=c(1,2),antigen_id=c(1,1)) 
antigen_map
```

    ## # A tibble: 2 × 2
    ##   exposure_id antigen_id
    ##         <dbl>      <dbl>
    ## 1           1          1
    ## 2           2          1

# 1.4 Force of Exposure and Exposure Model

Now, we specify the foe_pars argument which contains the force of
exposure for all exposure_IDs across all time steps. We also specify the
exposure model which will determine whether an individual is exposed to
a specific exposure event.

Since we did not specify different groups within demography, all
individuals will automatically be assigned group 1. Therefore, we only
need 1 row for dimension 1 in foe_pars. Groups can be used as an
indicator of location if the user wishes to specify a location specific
force of exposure. We specified the same value for all time steps within
foe_pars for simplicity but users will likely have varying numbers to
match real world settings.

``` r
## Create an empty array to store the force of infection for all exposure types
foe_pars <- array(0, dim=c(1,max(times),n_distinct(antigen_map$exposure_id)))
## Specify the force of exposure for exposure ID 1 which represents natural infection
foe_pars[,,1] <- 0.2
## Specify the force of exposure for exposure ID 2 which represents vaccination
foe_pars[,,2] <- 0.4

## Specify a simple exposure model which calculates the probability of exposure directly from the force of exposure at that time step. In this selected model, the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.
exposure_model<-exposure_model_simple_FOE
```

# 1.5 Immunity Model

Here, we specify the immunity model which will determine whether an
exposure event is successful or not. We will use a simple immunity
model(immunity_model_vacc_ifxn_simple) where successful exposure in only
conditional on the total number of previos exposure events. . With this
model, the probability of successful vaccination exposure depends on the
number of vaccines received prior to time t and age at time t while the
probability of successful infection is dependent on the number of
infections priot to time t.

``` r
## Specify immunity model within runserosim function below 
immunity_model<-immunity_model_vacc_ifxn_simple

## Specify 3 additional arguments needed for this immunity model 
## Specify which exposure IDs represent vaccination events 
vacc_exposures<-2
## Specify the age at which an individual is eligible for vaccination (9 months old); note non vaccine exposures are listed as NAs
vacc_age<-c(NA,9)
## Specify the maximum number of successful events an individual can receive for each exposure type (1 vaccine and 1 1 infection event)
max_events<-c(1,1)
```

# 1.6 Antibody Model and Model Parameters

Now, we specify the antibody model to be used within runserosim. We will
be using a monophasic boosting-waning model. This model assumes that for
each exposure there is a boost and boost waning parameter.The antibody
kinetics parameters are pre-loaded within a csv file. Users can edit the
csv file to specify their own parameters. All parameters needed for the
user specified antibody model, must be specified within the antibody
kinetics parameters tibble (theta). Lastly, we define the
draw_parameters function which determines how each individual’s antibody
kinetics parameters are simulated from the antibody kinetics parameters
tibble (theta). We will use a function which draws parameters directly
from theta for the antibody model with random effects. Parameters are
drawn randomly from a distribution with mean and standard deviation
specified within theta.

``` r
## Specify the antibody model 
antibody_model<-antibody_model_monophasic

## Bring in the antibody parameters needed for the antibody model
## Note that the titer-mediated protection parameters needed for the immunity model (Section 1.5), the titer-ceiling parameters needed for draw_parameters and the observation error parameter needed for the observation model (Section 1.7) are all defined here too.
## Also note that these are all arbitrary parameter values loosely informed by plausible values.
theta_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
theta <- read.csv(file = theta_path, header = TRUE)
theta
```

    ##   exposure_id antigen_id   name   mean     sd distribution
    ## 1           1          1  boost 4.0000 2.0000   log-normal
    ## 2           1          1   wane 0.0033 0.0005   log-normal
    ## 3          NA          1 obs_sd     NA 0.2500       normal
    ## 4           2          1  boost 2.0000 1.0000   log-normal
    ## 5           2          1   wane 0.0016 0.0005   log-normal

``` r
## Specify the draw_parameters function
draw_parameters<-draw_parameters_random_fx
```

# 1.7 Observation Model and observation_times

Now we specify how observed antibody titers are generated as a
probabilistic function of the true, latent antibody titer and when to
observe these titers. In this step, we specify the sampling design and
assay choice for their serological survey.We will take samples of all
individuals the end of the simulation (t=120).

Our chosen observation model observes the latent titer values given a
continuous assay with added noise. The added noise represents assay
variability and is implemented by sampling from a distribution with the
true latent antibody titer as the mean and the measurement error as the
standard deviation. The observation standard deviation and distribution
are defined within theta as the “obs_sd” parameter.

``` r
## Specify the observation model 
observation_model<-observation_model_continuous_noise

## Specify observation_times (serological survey sampling design) to observe antigen 1 across all individuals at the end of the simulation (t=120)
observation_times<- tibble(i=1:max(demography$i),t=120, ag=1)
```

# 1.8 Optional Arguments

There are no optional arguments needed for this simulation.

# 1.9 Run Simulation

This is the core simulation where all simulation settings, models and
parameters are specified within the main simulation function. Time to
run this step varies depending on the number of individuals and the
complexities of the specified models.

``` r
## Run the core simulation and save outputs in "res"
res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  antigen_map,
  theta,
  exposure_model,
  immunity_model,
  antibody_model,
  observation_model,
  draw_parameters,

  ## Other arguments needed
  max_events=max_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
)

## Note that models and arguments specified earlier in the code can be specified directly within this function.
```

# 1.10 Generate Plots

Now that the simulation is complete, let’s plot and examine the
simulation outputs.

``` r
## Plot antibody states and exposure histories for 10 individuals 
plot_subset_individuals_history(res$antibody_states, res$exposure_histories_long, subset=10, demography)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
## Plot exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)
```

![](README_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
## Plot exposure probabilities for all exposure types
plot_exposure_prob(res$exposure_probabilities_long)
```

![](README_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
## Plot antibody states for all individuals
plot_titers(res$antibody_states)
```

![](README_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

``` r
## Plot the serosurvey results 
plot_obs_titers_one_sample(res$observed_antibody_states)
```

![](README_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

``` r
## Note that the simulated kinetics parameters are also stored
head(res$kinetics_parameters)
```

    ## # A tibble: 6 × 7
    ##       i     t     e    ag name    value realized_value
    ##   <int> <dbl> <dbl> <int> <chr>   <dbl>          <dbl>
    ## 1     1   118     1     1 boost 4.39           4.39   
    ## 2     1   118     1     1 wane  0.00288        0.00288
    ## 3     2     7     1     1 boost 4.18           4.18   
    ## 4     2     7     1     1 wane  0.00463        0.00463
    ## 5     2    16     2     1 boost 1.72           1.72   
    ## 6     2    16     2     1 wane  0.00225        0.00225
