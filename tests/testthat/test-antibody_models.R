test_that("Check that antibody_model_monophasic function works", {
  
  ## Load in example data and necessary arguments
  tmp_pars <- list()
  tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,1,NULL, NULL, example_model_pars_numeric)
  
  ## Expect that the output of antibody_model_monophasic function is 4
  use_exposure_history <- xtabs(value ~ i + t + x, data=example_exposure_histories)
  use_biomarker_states <- xtabs(value ~ i + t + b, data=example_biomarker_states)
  use_exposure_history[1,1,1] <- 1
  expect_equal(antibody_model_monophasic(1,1,1,use_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric),4)
})

test_that("Check that antibody_model_monophasic and antibody_model_monophasic_cpp return the same values", {
  
  ## Load in example data and necessary arguments
  tmp_pars <- list()
  tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,1,NULL, NULL, example_model_pars_numeric)
  
  ## Expect that the output of antibody_model_monophasic function is 4
  use_exposure_history <- xtabs(value ~ i + t + x, data=example_exposure_histories)
  use_biomarker_states <- xtabs(value ~ i + t + b, data=example_biomarker_states)
  use_exposure_history[1,1,1] <- 1
  expect_equal(
    sapply(1:100, function(x) antibody_model_monophasic(1,x,1,use_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric)),
    sapply(1:100, function(x) antibody_model_monophasic_cpp(1,x,1,use_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric))
    )
})

test_that("Check that antibody_model_biphasic and antibody_model_biphasic_cpp return the same values", {
  
  ## Load in example data and necessary arguments
  tmp_pars <- list()
  tmp_model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,1,NULL, NULL, tmp_model_pars)
  
  ## Expect that the output of antibody_model_monophasic function is 4
  use_exposure_history <- xtabs(value ~ i + t + x, data=example_exposure_histories)
  use_biomarker_states <- xtabs(value ~ i + t + b, data=example_biomarker_states)
  use_exposure_history[1,1,1] <- 1
  expect_equal(
    sapply(1:100, function(x) antibody_model_biphasic(1,x,1,use_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric)),
    sapply(1:100, function(x) antibody_model_biphasic_cpp(1,x,1,use_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric))
  )
})


test_that("Check that antibody_model_biphasic function works", {
  
  ## Load in example data and necessary arguments 
  model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
 
  example_biomarker_map_numeric_tmp <- tibble(exposure_id = c(1,2,3,3),biomarker_id=c(1,2,1,2))
  ## Expect that the output of antibody_model_biphasic function is 0.075000975
  tmp_exposure_history <- array(0,dim=c(1,11,3))
  tmp_exposure_history[1,1,1] <- 1
  tmp_biomarker_states <- array(0, dim=c(1, 11, 3))
  
  tmp_pars <- list()
  tmp_pars[[1]] <- draw_parameters_fixed_fx_biomarker_dep(1,1,1,1,NULL, tmp_biomarker_states, model_pars)
  
  expect_equal(antibody_model_biphasic(1,1,1,tmp_exposure_history, use_biomarker_states, tmp_pars, example_biomarker_map_numeric_tmp),0.75)
})


test_that("Check that antibody_model_typhoid function works", {
 
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
