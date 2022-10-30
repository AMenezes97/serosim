library(microbenchmark)
microbenchmark(runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  antigen_map,
  theta,
  exposure_model=test,
  immunity_model,
  antibody_model,
  observation_model,
  draw_parameters,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  mod=mod,
  age_mod=age_mod,
  t_in_year=t_in_year
), 
runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  antigen_map,
  theta,
  exposure_model=test2,
  immunity_model,
  antibody_model,
  observation_model,
  draw_parameters,
  
  ## Other arguments needed
  boundary=boundary,
  max_vacc_events=max_vacc_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  mod=mod,
  age_mod=age_mod,
  t_in_year=t_in_year
), times=2)
