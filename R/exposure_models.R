#' Fixed probability of exposure
#' 
#' @description Generic wrapper function to allow the probability of exposure for each time point to be drawn directly from the foe_pars array.
#' @inheritParams exposure_model_simple_FOE
#' @param foe_pars A 3D array providing the *probability* of exposure for each exposure ID, group and time.
#' 
#' @return A probability of exposure from the foe_pars array
#' @export
#' @examples
#' times <- seq(0,365,by=1)
#' n_groups <- 1
#' n_exposures <- 1
#' foe_pars <- array(0.01, dim=c(n_groups,length(times),n_exposures))
#' exposure_model_fixed(1,1,1,1,foe_pars,NULL) 
exposure_model_fixed <- function(i, t, x, g, foe_pars, demography, ...){
   foe_pars[g, t, x]
}


#' Exposure Model Simple- Force of Exposure (FOE)
#' 
#' @description This is a simple exposure model where the probability of exposure depends on the force of exposure at that time for that group
#'
#' @param i individual
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
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements specified within the dem_mod. Within dem_mod, users can select which demographic elements affect the probability of exposure and by how much.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier (how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. Users can also add age modifier (how much each age affects probability of exposure). The column name will be "age" with the entry being individual's ages. 
#' @param t_in_year The number of time steps in a year; defaults to 1
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
exposure_model_dem_mod <- function(i, t, x, g, foe_pars, demography, dem_mod, t_in_year=1, ...){
  ## Find the force of exposure
  p <- foe_pars[g, t, x]
  
  ## Find which exposure events have modifiers
  exps<-unique(dem_mod$exposure_id)
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
    
    ## Pull out unique columns (demography variables) in dem_mod 
    cols<-unique(dem_mod$column)
    
    ## Pull individual's information within demography 
    individualnum<-i
    demography_tmp<-data.table(demography)
    demography_tmp<-demography_tmp[demography_tmp$i==individualnum & demography_tmp$times==t,]
    
    ## Convert dem_mod to data table
    mod2<-data.table(dem_mod)
    
    for (col in seq_along(cols)){ ## For each unique column entry in dem_mod 
      ## Pull the column name
      colname<-cols[col] 
      if(colname!="age"){
        ## Pull the individual's column entry within demography 
        entry<-demography_tmp[[colname]]
        ## Find the modifier within mod tibble 
        modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
        ## Multiply modifier by p 
        p <- p*modifier
      }
      if(colname=="age"){
        ## Calculate individual's current age
        ## Pull individual's birth time
        birth_time<-demography_tmp$birth
        curr_age<-floor((t-birth_time)/t_in_year)
        
        ## Find the age modifier within dem_mod tibble
        modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==curr_age]
        ## Multiply modifier by p 
        p <- p*modifier
      }
    }
    p_exp<-1-exp(-p)
    p_exp
  }
}

#' SIR Exposure Model
#' 
#' @description Finds the probability of exposure governed by an SIR model with specified parameters for each exposure type and group combination.
#' @inheritParams exposure_model_simple_FOE
#' @param foe_pars Data frame containing SIR model parameters for each group and exposure combination. Variable names: x (exposure ID), g (group ID), name (parameter name), value (parameter value). Parameters needed are: beta (transmission rate), gamma (recovery rate), I0 (per capita infected population seed size), R0 (per capita recovered population seed size) and t0 (seeding time). See example for format.
#' @param time_res Time steps to solve the ODEs. Set lower for higher accuracy.
#' @param ... Additional arguments.
#' @return Probability of exposure for the requested time step
#' @export
#' @examples 
#' times <- seq(0,365,by=1)
#' ## Create FOI (force of infection) from SIR model for one exposure type for one group
#' n_groups <- 1
#' n_exposures <- 1
#' ## Create parameters of the simple SIR model for one group and one exposure type
#' foe_pars <- data.frame(x=1,g=1,name=c("beta","gamma","I0","R0","t0"),values=c(0.2,1/7,1/10000,0,50))
#' ## Solve over all times as example
#' sir_prob <- exposure_model_sir(1, times, 1, 1, foe_pars)
#' plot_exposure_model(exposure_model_sir,seq(1,365,by=1),n_groups=1,n_exposures=1,foe_pars=foe_pars)
exposure_model_sir <- function(i, t, x, g, foe_pars, 
                               demography=NULL,
                               time_res=1,...){
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
    times <- seq(0,max(t),by=time_res)
    tmp_pars <- foe_pars[foe_pars$x == x & foe_pars$g == g,]
    pars <- tmp_pars$value
    names(pars) <- tmp_pars$name
    initial_states <- c(S=1-pars["I0"],I=pars["I0"],R=pars["R0"],inc=0)
    names(initial_states) <- c("S","I","R","inc")
    
    res <- deSolve::ode(y=initial_states, times=times, func=SIR_odes_with, parms= pars)
    diff(c(0,res[,"inc"]))[which(times %in% t)]
}

#' Simulate incidence curve using Gaussian Process
#'
#' @description Simulates an incidence curve (probability of infection per unit time) and associated parameters from a Gaussian Process model assuming that the covariance function (kernel) on time follows the squared exponential.
#' @param pars Vector of model parameters, must have entries 1) l: the lengthscale of the covariance function (length of the wiggles); 2) sigma: the output variance (average distance of function from mean); 3) scale_factor: parameter to scale the incidence curve to give a desired "overall probability" of infection (specifically, solves the model, transforms each entry to unit scale, and then scales the entire vector as scale_factor*p/(sum(p))); 4) tmax: the maximum time to solve the model over, generating a time vector starting at 0 and ending at tmax in increments of 1.
#' @return a list with two elements: 1) a vector giving the probability of infection at each timestep; 2) a vector giving the parameters used to solve the model
#' @examples
#' pars <- c("sigma"=1,"l"=100,"scale_factor"=1, "tmax"=365)
#' tmp <- simulate_gaussian_process(pars)
#' plot(tmp$incidence,type='l',ylim=c(0,0.01))
#' @export
simulate_gaussian_process <- function(pars){
    l <- unname(pars["l"])
    sigma <- unname(pars["sigma"])
    tmax <- unname(pars["tmax"])
    times <- seq(0,tmax,by=1)
    mus <- rep(0, length(times))
    
    ## Matrix of distances between each pair of times
    mat <- matrix(rep(times, each=length(times)),ncol=length(times))
    dist <- abs(apply(mat, 2, function(x) x-times))
    ## Generate covariance matrix
    K <- sigma^2 * exp(-(dist^2) / (2*l^2))
    ## Ensure positive definite
    diag(K) <- diag(K) + 1e-8
    L_K <- t(chol(K))
    
    ## Simulate
    eta <- rnorm(nrow(L_K),mean=0,sd=1)
    k <- (L_K %*% eta)[,1]
    ## Convert to probability
    ps <- 1/(1+exp(-k))
    
    scale_factor <- unname(pars["scale_factor"])
    ## Calculates incidence rate (probability of infection) for given times 
    ps <- (ps/sum(ps))*scale_factor
    
    list("incidence"=ps,"pars"=c("tmax"=tmax,"scale_factor"=scale_factor,"sigma"=sigma,"l"=l,"eta"=unname(eta)))
}


#' Gaussian process model
#' 
#' @description Generates an incidence curve (probability of infection per unit time) and associated parameters from a Gaussian Process model assuming that the covariance function (kernel) on time follows the squared exponential covariance function. It is recommended to use the outputs of \code{\link{simulate_gaussian_process}} as inputs to this function.
#' @inheritParams exposure_model_simple_FOE
#' @param foe_pars Data frame giving named parameters of the Gaussian process model for each unique group and exposure ID. Mandatory variable names are: 1) value (parameter value); 2) name (parameter name); 3) x (exposure ID); 4) g (group ID). Each group/exposure combination must have entries for the following parameters: l, sigma, scale_factor, tmax, eta. Note that eta gives the means of the Gaussian process at each time point, and thus there must be one entry for every element of seq(0,tmax,by=1).
#' 
#' @return Returns incidence for given time points t.
#' 
#' @examples 
#' pars_x1 <- c("sigma"=1,"l"=100,"scale_factor"=1, "tmax"=365*5)
#' pars_x2 <- c("sigma"=2,"l"=100,"scale_factor"=0.25, "tmax"=365*5)
#' tmp_x1 <- simulate_gaussian_process(pars_x1)
#' tmp_x2 <- simulate_gaussian_process(pars_x2)
#' foe_pars1 <- data.frame(name=names(tmp_x1$pars), value=unname(tmp_x1$pars),x=1,g=1)
#' foe_pars2 <- data.frame(name=names(tmp_x2$pars), value=unname(tmp_x2$pars),x=2,g=1)
#' foe_pars <- bind_rows(foe_pars1, foe_pars2)
#' exposure_model_gaussian_process(1, 365, 1, 1, foe_pars, NULL)
#' @export
exposure_model_gaussian_process <- function(i, t, x, g, foe_pars, demography, ...){
    ## Model parameters
    foe_pars <- foe_pars[foe_pars$x == x & foe_pars$g == g,]
    foe_pars_tmp <- foe_pars$value
    par_names <- foe_pars$name
    names(foe_pars_tmp) <- par_names
    
    times <- seq(0,foe_pars_tmp["tmax"],by=1)
    use_names <- c("eta",paste0("eta",1:length(times)))
    k <- foe_pars_tmp[which(par_names%in%use_names)]
    
    ## Create a matrix of times
    mat <- matrix(rep(times, each=length(times)),ncol=length(times))
    ## Subtract times from the columns of matrix mat and take the absolute value
    t_dist <- abs(apply(mat, 2, function(x) x-times))
    
    l <- foe_pars_tmp["l"]
    sigma <- foe_pars_tmp["sigma"]
    K <- sigma^2 * exp(-(dist^2) / (2*l^2))
    ## Add 0.01 to the diagonal of K to ensure positive definite
    diag(K) <- diag(K) + 1e-8
    
    ## Transpose the Choleski factorization of K
    L_K <- t(chol(K))
    
    ## Matrix multiplication
    k1 <- (L_K %*% k)
    ps <- 1/(1+exp(-k1))
    
    scale_factor <- unname(foe_pars_tmp["scale_factor"])

    ## Calculates incidence rate (probability of infection) for given times 
    prob_infection_tmp <- scale_factor*ps[,1]/sum(ps[,1])
    prob_infection_tmp[t]
}
