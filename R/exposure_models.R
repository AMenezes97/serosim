#' Exposure Model Simple- Force of Exposure (FOE)
#' 
#' @description This is a simple exposure model where the probability of exposure depends on the force of exposure at that time for that group
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
#' times <- seq(1,365,by=1)
#' ## Create fixed FOI (force of infection) for one exposure type for two groups
#' n_groups <- 2
#' n_exposures <- 1
#' foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
#' foe_pars[1,,] <- 0.01
#' foe_pars[2,,] <- 0.005
#' ## Solve the model for each time point for each group. Note that the i argument is not used, 
#' ## but is kept for compatibility with other functions.
#' foe <- matrix(NA, nrow=n_groups,ncol=length(times))
#' for(g in 1:nrow(foe)){
#'     foe[g,] <- unlist(sapply(times, function(t) exposure_model_simple_FOE(NULL, t, 1, g, foe_pars, NULL)))
#' }
exposure_model_simple_FOE <- function(i, t, x, g, foe_pars, demography, ...){
  p <- foe_pars[g, t, x]
  p_exp<-1-exp(-p)
  p_exp
}

#' Exposure Model Modified By Relevant Demographic Elements 
#'  
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements specified within the mod. Within mod, users can select which demographic elements affect the probability of exposure and by how much.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param ... 
#'  
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
#' times <- seq(1,365,by=1)
#' ## Create fixed FOI (force of infection) for two exposure types in one group
#' n_groups <- 1
#' n_exposures <- 2
#' n_times <- length(times)
#' n_indiv <- 2
#' foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
#' foe_pars[1,,1] <- 0.01
#' foe_pars[1,,2] <- 0.005
#' 
#' ## Create demography modifiers
#' ## Example with two individuals, one in low SES and one in high SES
#' demography <- tibble(i = rep(1:n_indiv, each=n_times), t=rep(times,2),SES=rep(c("low","high"),each=n_times))
#' 
#' ## Create example where for exposure ID 1, high SES gives 25% reduction in FOE.
#' ## high SES gives 50% reduction in FOE for exposure ID 2
#' dem_mod <- tibble(exposure_id=c(1,1,2,2),column=c("SES","SES","SES","SES"),
#'                  value=c("low","high","low","high"),modifier=c(1,0.75,1,0.5))
#' 
#' ## Solve the model for each time point for each group.
#' ## but is kept for compatibility with other functions.
#' foe <- array(NA, dim=c(n_indiv, n_times, n_exposures))
#' for(i in 1:n_indiv){
#'     for(x in 1:n_exposures){
#'        foe[i,,x] <- unlist(sapply(times, function(t) exposure_model_dem_mod(i, t, x, 1, foe_pars, demography, dem_mod)))
#'     }
#' }
exposure_model_dem_mod <- function(i, t, x, g, foe_pars, demography, dem_mod, ...){
  ## Find the force of exposure
  p <- foe_pars[g, t, x]
  
  ## Find which exposure events have modifiers
  exps<-dem_mod$exposure_id %>% unique()
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
      ## Pull out unique columns (demography variables) in dem_mod 
      cols<-dem_mod$column %>% unique()
      
      ## Pull individual's information within demography 
      demography_tmp<-data.table(demography)
      demography_tmp<-demography_tmp[demography_tmp$i==i & times==t,]
      
      for (col in seq_along(cols)){ ## For each unique column entry in dem_mod 
        ## Pull the column name
        colname<-cols[col] 
        ## Pull the individual's column entry within demography 
        entry<-demography_tmp[[colname]]
        ## Find the modifier within mod tibble 
        mod2<-data.table(dem_mod)
        modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
        ## Multiply modifier by p 
        p <- p*modifier
      }
  }
  p_exp<-1-exp(-p)
  p_exp
}


#' Exposure Model Modified By Relevant Demographic Elements and Age
#' 
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements (mod) and age (age_mod). This model combines exposure_model_V2 and exposure_model_V3
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param age_mod A tibble specifying the age modifier(how much each age affects probability of exposure); column names are age and modifier. All ages in simulation must have a modifier value within this tibble.
#' @param t_in_year The number of time steps in a year; defaults to 1 
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_dem_age_mod <- function(i, t, x, g, foe_pars, demography, dem_mod, age_mod, t_in_year=1, ...){
  ## Find the force of exposure 
  p <- foe_pars[g, t, x]
  
  ## Find which exposure events have modifiers
  exps<-dem_mod$exposure_id %>% unique()
  
  ## Pull individual's information within demography 
  demography_tmp<-data.table(demography)
  demography_tmp<-demography_tmp[demography_tmp$i==i & times==t,]
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
    
    ## Pull out unique columns in den_mod 
    cols<-dem_mod$column %>% unique()
    
    for (column in seq_along(cols)){ ## For each unique column entry in dem_mod 
      ## Pull the column name
      colname<-cols[column] 
      ## Pull the individual's column entry within demography 
      entry<-demography_tmp[[colname]]
      ## Find the modifier within mod tibble 
      mod2<-data.table(dem_mod)
      modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
      ## Multiply modifier by p
      p <- p*modifier
    }
  }
  
  ## Check if age modifiers apply to this particular exposure 
  exps_age<-age_mod$exposure_id %>% unique()
  
  ## If there are defined modifiers for that exposures then
  if(x  %in% exps_age){ 
    
    age_mod2<-data.table(age_mod)
    age_mod_tmp<-age_mod2[age_mod$exposure_id==x,]
    
    ## Calculate individual's current age
    birth_time<-demography_tmp$birth
    curr_age<-floor((t-birth_time)/t_in_year) ## Individual's birth time already gets pulled earlier on within the runserosim code so is it fine to just call it here?
    
    ## Find the age modifier within age_mod tibble 
    age_modifier<-age_mod_tmp$modifier[age_mod_tmp$age==curr_age]
    
    ## Multiply age modifier by p
    p <- p*age_modifier
  }
    p_exp<-1-exp(-p)
    p_exp
  
}
  
#' SIR Exposure Model
#' 
#' @inheritParams exposure_model_simple_FOE
#' @param foe_pars Data frame containing SIR model parameters for each group and exposure combination. Variable names: x (exposure ID), g (group ID), name (parameter name), value (parameter value). Parameters needed are: beta (transmission rate), gamma (recovery rate), I0 (per capita infected population seed size), R0 (per capita recovered population seed size) and t0 (seeding time).
#' @param time_res Time steps to solve the ODEs. Set lower for higher accuracy.
#' @param ...
#' @return Probability of exposure for the requested time step
#' @export
#' @examples 
exposure_model_sir <- function(i, t, x, g, foe_pars, demography=NULL,time_res=1,...){
    SIR_odes_with <- function(t, y, pars){
        with(as.list(c(y,pars)),{
            if(t > t0){
                dS <- -beta*S*I
                dI <- beta*S*I - gamma*I
                dR <- gamma*I
                dInc <- beta*S*I
            } else {
                dS <- 0
                dI <- 0
                dR <- 0
                dInc <- 0
            }
            
            return(list(c(dS, dI, dR, dInc)))
        })
    }
    times <- seq(0,t,by=time_res)
    tmp_pars <- foe_pars[foe_pars$x == x & foe_pars$g == g,]
    pars <- tmp_pars$value
    names(pars) <- tmp_pars$name
    initial_states <- c(S=1-pars["I0"],I=pars["I0"],R=pars["R0"],inc=0)
    names(initial_states) <- c("S","I","R","inc")
    
    res <- deSolve::ode(y=initial_states, times=times, func=SIR_odes_with, parms= pars)
    diff(c(0,res[,"inc"]))[which(times == t)]
}

