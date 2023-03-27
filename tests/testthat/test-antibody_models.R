test_that("Check that antibody_model_monophasic function works", {
  
  ## Load in example data and necessary arguments
  tmp_pars <- list()
  tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,1,NULL, NULL, example_model_pars_numeric)
  
  ## Expect that the output of antibody_model_monophasic function is 4
  expect_equal(antibody_model_monophasic(1,1,1,example_exposure_histories, example_biomarker_states, tmp_pars, example_biomarker_map_numeric),4)
})


test_that("Check that antibody_model_biphasic function works", {
  
  ## Load in example data and necessary arguments 
  model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  tmp_pars <- list()
  tmp_pars[[1]] <- draw_parameters_fixed_fx_biomarker_dep(1,1,1,1,NULL, NULL, model_pars)
  
  ## Expect that the output of antibody_model_biphasic function is 0.075000975
  expect_equal(antibody_model_biphasic(1,1,1,example_exposure_histories, example_biomarker_states, tmp_pars, example_biomarker_map_numeric),0.075000975)
})


test_that("Check that antibody_model_biphasic function works", {
 
   ## Load in example data and necessary arguments
    model_pars <- read.csv(system.file("extdata", "model_pars_typhoid.csv", package="serosim"))
    tmp_pars <- list()
    tmp_pars[[1]] <- draw_parameters_random_fx(1,1,1,1,NULL,NULL,model_pars)
    tmp_exposure_history <- array(0,dim=c(1,11,2))
    tmp_exposure_history[1,1,1] <- 1
    outcome<-antibody_model_typhoid(1,10, 1, tmp_exposure_history, NULL, tmp_pars,example_biomarker_map_numeric)

  ## Expect that the output of antibody_model_biphasic function is 0.075000975
  expect_equal(is.numeric(outcome),TRUE)
})
