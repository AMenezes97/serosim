test_that("Check that immunity_model_all_successful function works", {
  
  ## Load in example data and necessary arguments
  
  ## Expect that the output of immunity_model_all_successful function is 1
  expect_equal(immunity_model_all_successful(1,1,1,1),1)
})


test_that("Check that immunity_model_vacc_only function works if individual is eligible for vaccination", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 1))
  tmp_exposure_history[1,1,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_only function is 1
  expect_equal(immunity_model_vacc_only(1,8,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_vacc_events=2,vacc_age=5),1)
})


test_that("Check that immunity_model_vacc_only function works if indivdiual has already had the max number of vaccines", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 1))
  tmp_exposure_history[1,1,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_only function is 1
  expect_equal(immunity_model_vacc_only(1,8,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_vacc_events=1,vacc_age=5),0)
})


test_that("Check that immunity_model_vacc_only function works if individual is ineligible because they are too young", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 1))
  tmp_exposure_history[1,1,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_only function is 1
  expect_equal(immunity_model_vacc_only(1,3,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_vacc_events=2,vacc_age=5),0)
})


test_that("Check that immunity_model_vacc_ifxn_simple function works if individual is ineligible for vaccination exposure ", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 early exposures with exposure ID 1
  tmp_exposure_history[1,1:3,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_ifxn_simple function is 1
  expect_equal(immunity_model_vacc_ifxn_simple(1,8,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_events=c(3,5),vacc_exposures=c(1,2),vacc_age=c(1,1)),0)
})


test_that("Check that immunity_model_vacc_ifxn_simple function works if individual is eligible for vaccination exposure ", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 early exposures with exposure ID 1
  tmp_exposure_history[1,1:3,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)

  ## Expect that the output of immunity_model_vacc_ifxn_simple function is 1
  expect_equal(immunity_model_vacc_ifxn_simple(1,8,2,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_events=c(3,5),vacc_exposures=c(1,2),vacc_age=c(1,1)),1)
})


test_that("Check that immunity_model_vacc_ifxn_simple function works if individual is ineligible for natural infection exposure ", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 early exposures with exposure ID 1
  tmp_exposure_history[1,1:3,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_ifxn_simple function is 0
  expect_equal(immunity_model_vacc_ifxn_simple(1,8,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_events=c(3,5),vacc_exposures=2,vacc_age=c(1,1)),0)
})


test_that("Check that immunity_model_vacc_ifxn_simple function works if individual is eligible for natural infection exposure ", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 early exposures with exposure ID 1
  tmp_exposure_history[1,1:3,1] <- 1
  tmp_demography <- tibble(i=1, birth=1)
  
  ## Expect that the output of immunity_model_vacc_ifxn_simple function is 1
  expect_equal(immunity_model_vacc_ifxn_simple(1,8,1,tmp_exposure_history,NULL, tmp_demography,NULL, NULL,max_events=c(4,5),vacc_exposures=2,vacc_age=c(1,1)),1)
})


test_that("Check that immunity_model_ifxn_biomarker_prot function works", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 1))
  ## Toy example: individual has 1 prior exposure
  tmp_exposure_history[1,1,1] <- 1
  ## Set all biomarker states to 3 for sake of example
  tmp_biomarker_states <- array(0, dim=c(1,10,1))
  tmp_biomarker_states[1,,1] <- 3
  tmp_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  
  ## Expect that the output of immunity_model_ifxn_biomarker_prot function is 0.4273919
  expect_equal(immunity_model_ifxn_biomarker_prot(1,8,1,exposure_histories=tmp_exposure_history, biomarker_states=tmp_biomarker_states, demography=NULL, biomarker_map=example_biomarker_map_numeric, model_pars=tmp_pars),0.42739194)
})


test_that("Check that immunity_model_vacc_ifxn_biomarker_prot function works for eligible vaccination exposures", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 prior exposures to exposure ID 1, and none to exposure ID 2
  tmp_exposure_history[1,1:3,1] <- 1
  ## Set all biomarker states to 3 for sake of example
  tmp_biomarker_states <- array(0, dim=c(1,10,1))
  tmp_biomarker_states[1,,1] <- 3
  tmp_demography <- tibble(i=1, birth=1)
  tmp_pars <- reformat_biomarker_map(example_model_pars_biphasic)

  ## Expect that the output of immunity_model_vacc_ifxn_biomarker_prot function is 1
  expect_equal(immunity_model_vacc_ifxn_biomarker_prot(1,8,1,exposure_histories=tmp_exposure_history, biomarker_states=tmp_biomarker_states, demography=tmp_demography, 
                                                       biomarker_map=example_biomarker_map_numeric, model_pars=tmp_pars,max_events=4,vacc_exposures=1,vacc_age=1)
,1)
})


test_that("Check that immunity_model_vacc_ifxn_biomarker_prot function works for ineligible vaccination exposures", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 prior exposures to exposure ID 1, and none to exposure ID 2
  tmp_exposure_history[1,1:3,1] <- 1
  ## Set all biomarker states to 3 for sake of example
  tmp_biomarker_states <- array(0, dim=c(1,10,1))
  tmp_biomarker_states[1,,1] <- 3
  tmp_demography <- tibble(i=1, birth=1)
  tmp_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  
  ## Expect that the output of immunity_model_vacc_ifxn_biomarker_prot function is 0
  expect_equal(immunity_model_vacc_ifxn_biomarker_prot(1,8,1,exposure_histories=tmp_exposure_history, biomarker_states=tmp_biomarker_states, demography=tmp_demography, biomarker_map=example_biomarker_map_numeric, 
                                                       model_pars=tmp_pars,max_events=3,vacc_exposures=1,vacc_age=1),0)
})


test_that("Check that immunity_model_vacc_ifxn_biomarker_prot function works for eligible natural infection exposures", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 prior exposures to exposure ID 1, and none to exposure ID 2
  tmp_exposure_history[1,1:3,1] <- 1
  ## Set all biomarker states to 3 for sake of example
  tmp_biomarker_states <- array(0, dim=c(1,10,1))
  tmp_biomarker_states[1,,1] <- 3
  tmp_demography <- tibble(i=1, birth=1)
  tmp_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  
  ## Expect that the output of immunity_model_vacc_ifxn_biomarker_prot function is  0.9999039
  expect_equal(immunity_model_vacc_ifxn_biomarker_prot(1,8,2,exposure_histories=tmp_exposure_history, biomarker_states=tmp_biomarker_states, demography=tmp_demography, biomarker_map=example_biomarker_map_numeric, model_pars=tmp_pars,max_events=c(3,10),vacc_exposures=1,vacc_age=1), 0.9999039)
})


test_that("Check that immunity_model_vacc_ifxn_biomarker_prot function works for ineligible natural infection exposures", {
  
  ## Load in example data and necessary arguments
  tmp_exposure_history <- array(0, dim=c(1, 10, 2))
  ## Toy example: individual has 3 prior exposures to exposure ID 1, and none to exposure ID 2
  tmp_exposure_history[1,1:3,2] <- 1
  ## Set all biomarker states to 3 for sake of example
  tmp_biomarker_states <- array(0, dim=c(1,10,1))
  tmp_biomarker_states[1,,1] <- 3
  tmp_demography <- tibble(i=1, birth=1)
  tmp_pars <- reformat_biomarker_map(example_model_pars_biphasic)
  
  ## Expect that the output of immunity_model_vacc_ifxn_biomarker_prot function is 0
  expect_equal(immunity_model_vacc_ifxn_biomarker_prot(1,8,2,exposure_histories=tmp_exposure_history, biomarker_states=tmp_biomarker_states, demography=tmp_demography, 
                                                       biomarker_map=example_biomarker_map_numeric, model_pars=tmp_pars,max_events=c(3,3),vacc_exposures=1,vacc_age=1),0)
})


