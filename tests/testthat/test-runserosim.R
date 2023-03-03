test_that("Check that runserosim works for README example", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,120,by=1) 
  demography <- generate_pop_demography(N=10, times=times, prob_removal=0)
  biomarker_map_original <- tibble(exposure_id=c("ifxn","vacc"),biomarker_id=c("IgG","IgG"))
  biomarker_map <-reformat_biomarker_map(biomarker_map_original)
  foe_pars <- array(0, dim=c(1,max(times),dplyr::n_distinct(biomarker_map$exposure_id)))
  foe_pars[,,1] <- 0.01
  foe_pars[,,2] <- 0.1
  model_pars_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
  model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
  model_pars<-reformat_biomarker_map(model_pars_original)
  
  expect_message(res<- runserosim(
    simulation_settings=list("t_start"=1,"t_end"=max(times)),
    demography,
    observation_times=tibble(i=1:max(demography$i),t=120, b=1),
    foe_pars,
    biomarker_map,
    model_pars,
    exposure_model=exposure_model_simple_FOE,
    immunity_model=immunity_model_vacc_ifxn_simple,
    antibody_model=antibody_model_monophasic,
    observation_model=observation_model_continuous_noise,
    draw_parameters=draw_parameters_random_fx,
    
    ## Other arguments needed
    max_events=c(1,1),
    vacc_exposures=2,
    vacc_age=c(NA,9),
    sensitivity=0.85,
    specificity=0.9
  ),"")
  
})


test_that("Check that runserosim works for case study 1 example", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,120,by=1) 
  demography <- generate_pop_demography(10, times, age_min=0, prob_removal=0)
  biomarker_map_original <- tibble(exposure_id=c("measles_ifxn","measles_vacc"),biomarker_id=c("measles_IgG","measles_IgG"))
  biomarker_map <-reformat_biomarker_map(biomarker_map_original)
  foe_pars <- array(0, dim=c(1,max(times),dplyr::n_distinct(biomarker_map$exposure_id)))
  foe_pars[,,1] <- 0.02
  foe_pars[,,2] <- 0.04
  model_pars_path <- system.file("extdata", "model_pars_cs1.csv", package = "serosim")
  model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
  model_pars<-reformat_biomarker_map(model_pars_original)
  obs1 <- tibble(i=1:10,t=60, b=1)
  obs2 <- tibble(i=1:10,t=120, b=1)
  observation_times<-rbind(obs1,obs2)
  
  expect_message(res<- runserosim(
    simulation_settings=list("t_start"=1,"t_end"=max(times)),
    demography,
    observation_times=observation_times,
    foe_pars,
    biomarker_map,
    model_pars,
    exposure_model=exposure_model_simple_FOE,
    immunity_model=immunity_model_vacc_ifxn_biomarker_prot,
    antibody_model=antibody_model_biphasic,
    observation_model=observation_model_continuous_bounded_noise,
    draw_parameters=draw_parameters_random_fx_biomarker_dep,
    
    ## Other arguments needed
    max_events=c(1,1),
    vacc_exposures=2,
    vacc_age=c(NA,9),
    sensitivity=0.996,
    specificity=1,
    bounds=tibble(biomarker_id=c(1,1),name=c("lower_bound","upper_bound"),value=c(100,Inf)),
    VERBOSE=5
  ),"")
  
})


test_that("Check that runserosim works for case study 2 example", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,120,by=1) 
  aux <- list("NS"=list("name"="NS","options"=c("low","medium","high"), "distribution"=c(0.3,0.3,0.4)),
              "Group"=list("name"="group","options"=c(1, 2), "distribution"=c(0.5,0.5)))
  demography <- generate_pop_demography(10, times, age_min=0, removal_min=1, removal_max=120, prob_removal=0.2, aux=aux)
  biomarker_map_original <- tibble(exposure_id=c("DP_ifxn","PT_ifxn","vacc","vacc"),biomarker_id=c("DP_antibody","PT_antibody","DP_antibody","PT_antibody"))
  biomarker_map <-reformat_biomarker_map(biomarker_map_original)
  foe_pars <- array(0, dim=c(dplyr::n_distinct(demography$group),max(times), dplyr::n_distinct(biomarker_map$exposure_id)))
  foe_pars[1,,1] <- 0.04
  foe_pars[2,,1] <- 0.03
  foe_pars[1,,2] <- 0.02 
  foe_pars[2,,2] <- 0.01 
  foe_pars[1,,3] <- 0.02 
  foe_pars[2,,3] <- 0.03
  age_mod_1<-tibble(exposure_id=rep(1,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))
  age_mod_2<-tibble(exposure_id=rep(2,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,1,1,1,1,1,1,1,1))
  age_mod_3<-tibble(exposure_id=rep(3,11), column=rep("age",times=11),  value=0:10, modifier=c(3,1,1,1,1,1,1,1,1,1,1))
  age_mod<-rbind(age_mod_1,age_mod_2,age_mod_3)
  mod<-tibble(exposure_id=c(1,1,1,2,2,2,3,3,3), column=rep("NS",times=9), value=rep(c("low","medium", "high"),3), modifier=c(2,1.5,1,2,1.5,1,1,2,3))
  dem_mod<-rbind(age_mod,mod)
  model_pars_path <- system.file("extdata", "model_pars_cs2.csv", package = "serosim")
  model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
  model_pars<-reformat_biomarker_map(model_pars_original)
  obs1 <- tibble(i=1:10,t=120, b=1)
  obs2 <- tibble(i=1:10,t=120, b=2)
  observation_times<-rbind(obs1,obs2)
  
  expect_message(res<- runserosim(
    simulation_settings=list("t_start"=1,"t_end"=max(times)),
    demography,
    observation_times=observation_times,
    foe_pars,
    biomarker_map,
    model_pars,
    exposure_model=exposure_model_dem_mod,
    immunity_model=immunity_model_vacc_ifxn_biomarker_prot,
    antibody_model=antibody_model_biphasic,
    observation_model=observation_model_continuous_bounded_noise,
    draw_parameters=draw_parameters_random_fx_biomarker_dep,
    
    ## Other arguments needed
    max_events=c(Inf,Inf,3),
    vacc_exposures=3,
    vacc_age=c(NA,NA,2),
    bounds=tibble(biomarker_id=c(1,1,2,2),name=rep(c("lower_bound","upper_bound"),2),value=c(0.01,2,5,200)),
    dem_mod=dem_mod,
    t_in_year=12
  ),"")
  
})
