# library(tictoc)
# 
# ## Run Case Study 1 100 times 
# run_time<-data.frame(
#   run=1:100,
#   time=NA
# )
#   
# for (runs in 1:nrow(run_time)){
#   tic()
#   res<- runserosim(
#     simulation_settings,
#     demography,
#     observation_times,
#     foe_pars,
#     biomarker_map,
#     model_pars,
#     exposure_model,
#     immunity_model,
#     antibody_model,
#     observation_model,
#     draw_parameters,
#     
#     ## Specify other arguments needed
#     VERBOSE=VERBOSE,
#     bounds=bounds,
#     max_events=max_events,
#     vacc_exposures=vacc_exposures,
#     vacc_age=vacc_age,
#     sensitivity=sensitivity,
#     specificity=specificity
#   )
#   time<-toc()
#   time$callback_msg
#   run_time$time[runs]<- time$callback_msg
# }
# 
# run_time_100<-run_time
# write.csv(run_time_100,"/Users/arthurmenezes/downloads//run_time_100.csv", row.names = FALSE)
#                                                                                                    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



# library(tictoc)
# 
# ## Run 100 times 
# run_time<-data.frame(
#   run=1:100,
#   time=NA
# )
# 
# for (runs in 1:nrow(run_time)){
#   tic()
#   res<- runserosim(
#     simulation_settings,
#     demography,
#     observation_times,
#     foe_pars,
#     biomarker_map,
#     model_pars,
#     exposure_model,
#     immunity_model,
#     antibody_model,
#     observation_model,
#     draw_parameters,
#     
#     ## Other arguments needed
#     max_events=max_events,
#     vacc_exposures=vacc_exposures,
#     vacc_age=vacc_age,
#     sensitivity=sensitivity,
#     specificity=specificity
#   )
#   time<-toc()
#   time$callback_msg
#   run_time$time[runs]<- time$callback_msg
# }
# 
# README_run_time_100<-run_time
# write.csv(README_run_time_100,"/Users/arthurmenezes/downloads//README_run_time_100.csv", row.names = FALSE)
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



library(tictoc)

## Run 100 times 
run_time<-data.frame(
  run=1:100,
  time=NA
)

for (runs in 1:nrow(run_time)){
  tic()
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
    dem_mod=dem_mod,
    t_in_year=t_in_year
  )
  time<-toc()
  time$callback_msg
  run_time$time[runs]<- time$callback_msg
}

cs2_run_time_250<-run_time
write.csv(cs2_run_time_250,"/Users/arthurmenezes/downloads//cs2_run_time_250.csv", row.names = FALSE)

























