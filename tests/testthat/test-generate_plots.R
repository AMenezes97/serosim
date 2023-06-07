test_that("Check that plot_biomarker_mediated_protection function works", {
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_biomarker_mediated_protection(0:10,5,0.9), regexp=NA)
})


test_that("Check that plot_biomarker_dependent_boosting function works", {
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_biomarker_dependent_boosting(0,10,1,4,.14), regexp=NA)
})


test_that("Check that plot_exposure_prob function works", {
  
  ## Load in example data 
  data("example_exposure_probabilities")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_exposure_prob(example_exposure_probabilities), regexp=NA)
})


test_that("Check that plot_exposure_force function works", {
  
  ## Load in example data 
  data("example_exposure_force")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_exposure_force(example_exposure_force), regexp=NA)
})


test_that("Check that plot_immune_histories function works", {
  
  ## Load in example data 
  data("example_immune_histories")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_immune_histories(example_immune_histories), regexp=NA)
})


test_that("Check that plot_biomarker_quantity function works", {
  
  ## Load in example data 
  data("example_biomarker_states")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_biomarker_quantity(example_biomarker_states), regexp=NA)
})


test_that("Check that plot_obs_biomarkers_one_sample function works", {
  
  ## Load in example data 
  data("example_observed_biomarker_states")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_obs_biomarkers_one_sample(example_observed_biomarker_states), regexp=NA)
})


test_that("Check that plot_obs_biomarkers_paired_sample function works", {
  
  ## Load in example data 
  data("example_biomarker_states")
  
  example_biomarker_states$observed <- example_biomarker_states$value
  example_biomarker_states_subset <- example_biomarker_states %>% filter(t %in% c(1,120))
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_obs_biomarkers_paired_sample(example_biomarker_states_subset), regexp=NA)
})


test_that("Check that plot_subset_individuals_history function works", {
  
  ## Load in example data 
  data("example_biomarker_states")
  data("example_immune_histories")
  data("example_demography")
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_subset_individuals_history(example_biomarker_states,example_immune_histories,3,example_demography), regexp=NA)
})


test_that("Check that plot_antibody_model function works with fixed effects (draw_parameters_fixed_fx)", {
  
  ## Load in example data and necessary arguments
  model_pars <-  read.csv(system.file("extdata", "model_pars_test_1.csv", package="serosim")) %>% drop_na()
  biomarker_map <- tibble(exposure_id=c(1,1,2),biomarker_id=c(1,2,1))
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_antibody_model(antibody_model_biphasic, 50, model_pars=model_pars,draw_parameters_fn = draw_parameters_fixed_fx, biomarker_map=biomarker_map), regexp=NA)
})


test_that("Check that plot_antibody_model function works with random effects (draw_parameters_random_fx)", {
  
  ## Load in example data and necessary arguments
  model_pars <- read.csv(system.file("extdata", "model_pars_test_1.csv", package="serosim")) %>% drop_na()
  biomarker_map <- tibble(exposure_id=c(1,1,2),biomarker_id=c(1,2,1))
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_antibody_model(antibody_model_biphasic, 50, model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx, biomarker_map=biomarker_map), regexp=NA)
})


test_that("Check that plot_exposure_model function works with case study 1 parameters", {
  
  ## Load in example data and necessary arguments
  ## Specify the number of time periods to simulate 
  times <- seq(1,120,by=1) 
  ## Create an empty array to store the force of exposure for all exposure types
  foe_pars <- array(0, dim=c(1,max(times),2))
  ## Specify the force of exposure for exposure ID 1 
  foe_pars[,,1] <- 0.02
  ## Specify the force of exposure for exposure ID 2 
  foe_pars[,,2] <- 0.04
  
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_exposure_model(exposure_mode=exposure_model_simple_FOE, times=times, n_groups = 1, n_exposures = 2, foe_pars=foe_pars), regexp=NA)
})


test_that("Check that plot_exposure_model function works with case study 2 parameters", {
  
  ## Load in example data and necessary arguments
  ## Specify the number of time periods to simulate 
  times <- seq(1,120,by=1) 
  
  ## Generate the population demography tibble
  aux <- list("NS"=list("name"="NS","options"=c("low","medium","high"), "distribution"=c(0.3,0.3,0.4)),
              "Group"=list("name"="group","options"=c(1, 2), "distribution"=c(0.5,0.5)))
  demography <- generate_pop_demography(100, times, age_min=0, removal_min=1, removal_max=120, prob_removal=0.2, aux=aux)
  ## Create an empty array to store the force of exposure for all exposure types
  foe_pars <- array(0, dim=c(2,max(times), 3))
  
  ## Specify the force of exposure for exposure ID 1 
  foe_pars[1,,1] <- 0.04 
  foe_pars[2,,1] <- 0.03 
  
  ## Specify the force of exposure for exposure ID 2 
  foe_pars[1,,2] <- 0.02 
  foe_pars[2,,2] <- 0.01 
  
  ## Specify the force of exposure for exposure ID 3 
  foe_pars[1,,3] <- 0.02
  foe_pars[2,,3] <- 0.03
  
  ## First, we will specify age modifiers 
  age_mod_1<-tibble(exposure_id=rep(1,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))
  age_mod_2<-tibble(exposure_id=rep(2,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,1,1,1,1,1,1,1,1))
  age_mod_3<-tibble(exposure_id=rep(3,11), column=rep("age",times=11),  value=0:10, modifier=c(3,1,1,1,1,1,1,1,1,1,1))
  age_mod<-rbind(age_mod_1,age_mod_2,age_mod_3)

  
  ## Specify additional demography exposure modifiers and combine them with the previous ones
  mod<-tibble(exposure_id=c(1,1,1,2,2,2,3,3,3), column=rep("NS",times=9), value=rep(c("low","medium", "high"),3), modifier=c(2,1.5,1,2,1.5,1,1,2,3))
  ## Combine both age modifiers and additional modifiers
  dem_mod<-rbind(age_mod,mod)
  
  ## Specify the number of time steps within a year which will be used to calculate an individual's age. 
  t_in_year=12
  
  
  ## Expect that plotting function works properly with no errors
  expect_message(plot_exposure_model(indivs=1:5,exposure_model=exposure_model_dem_mod, times=times, n_groups = 2, n_exposures = 3, foe_pars=foe_pars, demography=demography, dem_mod=dem_mod,t_in_year=t_in_year), regexp=NA)
})


test_that("Check that plot_exposure_model function works with first example listed in function documentation file", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,120,1)
  n_groups <- 1
  n_exposures <- 2
  foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
  foe_pars[1,,1] <- 0.01
  foe_pars[1,,2] <- 0.005 
  aux <- list("SES"=list("name"="SES","options"=c("low","high"), "distribution"=c(0.5,0.5)))
  demography <- generate_pop_demography(5, times, age_min=0, removal_min=1, removal_max=120, prob_removal=0.2, aux=aux)
  dem_mod <- tibble(exposure_id=c(1,1,2,2),column=c("SES","SES","SES","SES"),
                    value=c("low","high","low","high"),modifier=c(1,0.75,1,0.5))
  
  ## Expect that plotting function works properly with no errors
  expect_message( plot_exposure_model(indivs=1:4, exposure_model=exposure_model_dem_mod, times=times,1,2,foe_pars=foe_pars,demography = demography,dem_mod=dem_mod), regexp=NA)
})


test_that("Check that plot_exposure_model function works with second example (SIR model) listed in function documentation file", {
  
  ## Load in example data and necessary arguments
  ## SIR model with two groups and two exposure types
  foe_pars <- bind_rows(
                        tibble(x=1,g=1,name=c("beta","gamma","I0","R0","t0"),value=c(0.3,0.2,0.00001,0,0)),
                        tibble(x=2,g=1,name=c("beta","gamma","I0","R0","t0"),value=c(0.35,0.25,0.00001,0,200)),
                        tibble(x=1,g=2,name=c("beta","gamma","I0","R0","t0"),value=c(0.5,0.2,0.00005,0,0)),
                        tibble(x=2,g=2,name=c("beta","gamma","I0","R0","t0"),value=c(0.27,0.2,0.00001,0,50))
                        )

  ## Expect that plotting function works properly with no errors
  expect_message(plot_exposure_model(exposure_model=exposure_model_sir, times=seq(1,365,by=1),n_groups = 2,n_exposures = 2,foe_pars=foe_pars), regexp=NA)
})

