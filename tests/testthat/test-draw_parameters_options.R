test_that("Check that draw_parameters_fixed_fx function works", {
  
  ## Load in example data and necessary arguments
  draw_pars_tmp<-draw_parameters_fixed_fx(1,1,1,example_demography, example_biomarker_states, example_model_pars_numeric)
 
  ## Expect that the output of draw_parameters_fixed_fx function is 4
  expect_equal(draw_pars_tmp$realized_value[1], 4)
})


test_that("Check that draw_parameters_random_fx function works", {
  
  ## Load in example data and necessary arguments
 
  
  ## Expect that draw_parameters_random_fx function runs with no errors
  expect_message(draw_parameters_random_fx(1,1,1,example_demography, example_biomarker_states, example_model_pars_numeric), regexp=NA)
})


test_that("Check that draw_parameters_fixed_fx_biomarker_dep function works", {
  
  ## Load in example data and necessary arguments
  model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  draw_parameters_tmp<-draw_parameters_fixed_fx_biomarker_dep(2,100,1,example_demography, example_biomarker_states_wide, model_pars)
  
  ## Expect that the output of draw_parameters_fixed_fx_biomarker_dep function is 0.5
  expect_equal(draw_parameters_tmp$value[1], 0.5)
})


test_that("Check that draw_parameters_random_fx_biomarker_dep function works", {
  
  ## Load in example data and necessary arguments
  model_pars <- reformat_biomarker_map(example_model_pars_biphasic)

  ## Expect that draw_parameters_random_fx_biomarker_dep function runs with no errors
  expect_message( draw_parameters_random_fx_biomarker_dep(2,100,1,example_demography, example_biomarker_states_wide, model_pars), regexp=NA)
})
