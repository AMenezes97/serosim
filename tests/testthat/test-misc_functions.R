test_that("Check that generate_pop_demography function works with a default case", {

  ## Load in example data and necessary arguments
  dem_tmp<-generate_pop_demography(10, times=1:120, age_min=0, removal_min=1, removal_max=120, prob_removal=0.3)
  
  ## Expect that generate_pop_demography function works properly 
  expect_equal(max(dem_tmp$times), 120)
})


test_that("Check that generate_pop_demography function works when specifying birth_times", {
  
  ## Load in example data and necessary arguments
  birth_times <- rpois(100, 5)
  dem_tmp<-generate_pop_demography(N=100, times=1:120, birth_times=birth_times, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3)
  
  ## Expect that generate_pop_demography function works properly 
  expect_equal(max(dem_tmp$times), 120)
})


test_that("Check that generate_pop_demography function works when specifying additional demography elements", {
  
  ## Load in example data and necessary arguments
  aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "proportion"=c(0.5,0.5)),
              "Group"=list("name"="group","options"=c(1, 2, 3, 4), "proportion"=c(0.25,0.25,0.25,0.25)) )
  dem_tmp<-generate_pop_demography(10, 1:120, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3, aux=aux)
  col_names<-colnames(dem_tmp)
  
  ## Expect that generate_pop_demography function works properly 
  expect_equal(col_names[5], "sex")
})


test_that("Check that simulate_birth_times function works", {

  ## Load in example data and necessary arguments

  ## Expect that simulate_birth_times function works properly with no errors
  expect_message(simulate_birth_times(500, 1:100, age_min=9) , regexp = NA)
})


test_that("Check that simulate_removal_times  function works", {

  ## Load in example data and necessary arguments
  birth_times<-simulate_birth_times(500, 1:100, age_min=9) 

  ## Expect that simulate_removal_times  function works properly with no errors
  expect_message(simulate_removal_times(500,1:100,birth_times, removal_min=10,removal_max=99, prob_removal=0.4), regexp = NA)
})


test_that("Check that update function works", {
  
  ## Load in example data and necessary arguments
  
  ## Expect that update function works properly with no errors
  expect_message(update(10, 100), regexp = "")
})


test_that("Check that reformat_biomarker_map function works to convert characters to numeric", {
  
  ## Load in example data and necessary arguments
  biomarker_map <- tibble(exposure_id=c("infection","vaccination"),biomarker_id=c("IgG","IgG"))
  biomarker_map <- reformat_biomarker_map(biomarker_map)
  
  ## Expect that reformat_biomarker_map function works properly and produces numeric outputs
  expect_equal(is.numeric(biomarker_map$exposure_id),TRUE)
})


test_that("Check that reformat_biomarker_map function works to convert numbers to numeric", {
  
  ## Load in example data and necessary arguments
  biomarker_map <- tibble(exposure_id=c(1,2),biomarker_id=c(1,1))
  biomarker_map <- reformat_biomarker_map(biomarker_map, exposure_key=c("infection","vaccination"),biomarker_key=c(1))

  ## Expect that reformat_biomarker_map function works properly with no errors
  expect_equal(is.numeric(biomarker_map$exposure_id),FALSE)
})



test_that("Check that precomputation_checks function works", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,100,by=1)
  N <- 100
  n_exposure_ids <- 2
  n_groups <- 2
  demography <- generate_pop_demography(N=N, times=times,prob_removal=0,aux=list(Group=list(name="group",options=c(1,2),proportion=c(0.5,0.5))))
  foe_pars <- array(runif(n_exposure_ids*length(times)*n_groups),dim=c(n_groups,length(times),n_exposure_ids))
  res <- precomputation_checks(N, times=times, exposure_ids=1:2,groups=1:2,exposure_model_simple_FOE,foe_pars=foe_pars, demography=demography,VERBOSE=10,check_correct=TRUE)

  ## Expect that precomputation_checks function works properly with no errors and returns $flag == TRUE
  expect_equal(res$flag, TRUE)
})


test_that("Check that convert_indices_matrix_to_vector function works", {
  
  ## Load in example data and necessary arguments
  x <- matrix(runif(10*5),nrow=10,ncol=5)
  a<-x[3,3]
  b<-x[convert_indices_matrix_to_vector(3,3,10)]
  
  ## Expect that convert_indices_matrix_to_vector function works properly and produces equal values 
  expect_equal(a, b)
})


test_that("Check that normal_to_lognormal_mean function works", {
  
  ## Load in example data and necessary arguments
  
  ## Expect that normal_to_lognormal_mean function output equals 3.8924126
  expect_equal(normal_to_lognormal_mean(50,10), 3.8924126)
})



test_that("Check that normal_to_lognormal_sd function works", {

  ## Load in example data and necessary arguments

  ## Expect that normal_to_lognormal_sd function output equals 0.1980422
  expect_equal(normal_to_lognormal_sd(50,10), 0.1980422)
})
