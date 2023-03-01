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


