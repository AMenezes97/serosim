test_that("Check that exposure_model_fixed function works", {
 
  ## Load in example data and necessary arguments
  times <- seq(0,365,by=1)
  n_groups <- 1
  n_exposures <- 1
  foe_pars <- array(0.01, dim=c(n_groups,length(times),n_exposures))
  
  ## Expect that exposure probability is equal to 0.01
  expect_equal(exposure_model_fixed(1,1,1,1,foe_pars,NULL), 0.01)
})


test_that("Check that exposure_model_simple_FOE function works", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,365,by=1)
  n_groups <- 2
  n_exposures <- 1
  foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
  foe_pars[1,,] <- 0.01
  foe_pars[2,,] <- 0.005
  
  
  ## Expect that exposure probability is equal to 0.0099501663
  expect_equal(exposure_model_simple_FOE(1, 1, 1, 1, foe_pars, demography=NULL), 0.0099501663)
})


test_that("Check that exposure_model_dem_mod function works", {
  
  ## Load in example data and necessary arguments
  times <- seq(1,10,by=1)
  n_groups <- 1
  n_exposures <- 2
  n_times <- length(times)
  n_indiv <- 2
  foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
  foe_pars[1,,1] <- 0.01
  foe_pars[1,,2] <- 0.005

  ## Create demography modifiers
  ## Example with two individuals, one in low SES -and one in high SES
  demography <- tibble(i = rep(1:n_indiv, each=n_times), times=rep(times,2),SES=rep(c("low","high"),each=n_times))

  ## Create example where for exposure ID 1, high SES gives 25% reduction in FOE.
  ## high SES gives 50% reduction in FOE for exposure ID 2
  dem_mod <- tibble(exposure_id=c(1,1,2,2),column=c("SES","SES","SES","SES"),
                   value=c("low","high","low","high"),modifier=c(1,0.75,1,0.5))
  
  ## Expect that exposure probability is equal to 0.0074719452
  expect_equal(exposure_model_dem_mod(2, 1, 1, 1, foe_pars, demography=demography,dem_mod=dem_mod), 0.0074719452)
})


test_that("Check that exposure_model_dem_mod function works with age modifier", {
  
  ## Load in example data and necessary arguments
  ## Specify the number of time periods to simulate 
  times <- seq(1,10,by=1) 
  n_times <- length(times)
  n_indiv <- 2
  n_groups <- 1
  n_exposures <- 2
  
  ## Example with two individuals, one in low SES -and one in high SES
  demography <- tibble(i = rep(1:n_indiv, each=n_times), times=rep(times,2),birth=rep(1,20),SES=rep(c("low","high"),each=n_times))
  
  ## Create an empty array to store the force of exposure for all exposure types
  foe_pars <- array(0, dim=c(1,max(times), 2))
  
  ## Specify the force of exposure for exposure ID 1 
  foe_pars[1,,1] <- 0.04 
  
  ## Specify the force of exposure for exposure ID 2 
  foe_pars[1,,2] <- 0.02

  ## Specify age modifiers 
  age_mod_1<-tibble(exposure_id=rep(1,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,2,1,1,1,1,1,1,1))
  age_mod_2<-tibble(exposure_id=rep(2,11), column=rep("age",times=11), value=0:10, modifier=c(2,2,2,1,1,1,1,1,1,1,1))
  age_mod<-rbind(age_mod_1,age_mod_2)

  ## Specify additional demography exposure modifiers and combine them with the previous ones
  mod<-tibble(exposure_id=c(1,1,2,2), column=rep("SES",times=4), value=rep(c("low", "high"),2), modifier=c(2,1.5,1,2))
  ## Combine both age modifiers and additional modifiers
  dem_mod<-rbind(age_mod,mod)
  
  ## Specify the number of time steps within a year which will be used to calculate an individual's age. 
  t_in_year=1
  
  ## Expect that exposure probability is equal to 0.14785621
  expect_equal(exposure_model_dem_mod(1, 2, 1, 1, foe_pars, demography=demography,dem_mod=dem_mod, t_in_year=t_in_year), 0.14785621)
})


test_that("Check that exposure_model_sir function works", {
  
  ## Load in example data and necessary arguments
  times <- seq(0,365,by=1)
  n_groups <- 1
  n_exposures <- 1
  foe_pars <- data.frame(x=1,g=1,name=c("beta","gamma","I0","R0","t0"),values=c(0.2,1/7,1/10000,0,50))
  
  ## Expect that the exposure model function works properly with no errors
  expect_message(exposure_model_sir(1, times, 1, 1, foe_pars), regexp=NA)
})


test_that("Check that simulate_gaussian_process function works", {
  
  ## Load in example data and necessary arguments
  pars <- c("sigma"=1,"l"=100,"scale_factor"=1, "tmax"=365)
  
  ## Expect that the simulate_gaussian_process function works properly with no errors
  expect_message(simulate_gaussian_process(pars), regexp=NA)
})


test_that("Check that exposure_model_gaussian_process function works", {
  
  ## Load in example data and necessary arguments
  pars_x1 <- c("sigma"=1,"l"=100,"scale_factor"=1, "tmax"=365*5)
  pars_x2 <- c("sigma"=2,"l"=100,"scale_factor"=0.25, "tmax"=365*5)
  tmp_x1 <- simulate_gaussian_process(pars_x1)
  tmp_x2 <- simulate_gaussian_process(pars_x2)
  foe_pars1 <- data.frame(name=names(tmp_x1$pars), value=unname(tmp_x1$pars),x=1,g=1)
  foe_pars2 <- data.frame(name=names(tmp_x2$pars), value=unname(tmp_x2$pars),x=2,g=1)
  foe_pars <- bind_rows(foe_pars1, foe_pars2)
  
  ## Expect that the exposure_model_gaussian_process function works properly with no errors
  expect_message(exposure_model_gaussian_process(1, 365, 1, 1, foe_pars, NULL), regexp=NA)
})


test_that("Check that exposure_model_indiv_fixed function works", {
  
  ## Load in example data and necessary arguments
  times <- seq(0,365,by=1)
  n_indivs <- 100
  n_exposures <- 1
  foe_pars <- array(0.01, dim=c(n_indivs,length(times),n_exposures))
  
  ## Expect that the exposure model function works properly with no errors
  expect_equal(exposure_model_indiv_fixed(1,1,1,1,foe_pars,NULL), 0.01)
})


