#' Draw Parameters Fixed Effects
#'  
#' @description This function draws all parameters listed in theta for successful exposures with fixed effects.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param ag antigen
#' @param demography Demography information 
#' @param antibody_states An array of true antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_fixed_fx <- function(i, t, e, ag, demography, antibody_states, theta, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    pars[par] <- theta_tmp$mean[par]
    par_names[par] <- theta_tmp$name[par]
  }
  # all_pars <- tibble(i=i, t=t, e=e, ag=ag, name=par_names, value=pars) ## This line doesn't work
  all_pars <- tibble(i=rep(i,nrow(theta_tmp)),t=rep(t,nrow(theta_tmp)), e=rep(e,nrow(theta_tmp)), ag=theta_tmp$antigen_id, name=par_names, value=pars, realized_value=pars) 
  return(all_pars)
}
#' Draw Parameters Random Effects
#'  
#' @description This function draws all parameters listed in theta for successful exposures with random effects.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param ag antigen
#' @param demography Demography information 
#' @param antibody_states An array of true antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters  
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_random_fx_boost_wane <- function(i, t, e, ag, demography, antibody_states, theta, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    if(theta_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
      ## Create functions to convert normal distributions to log-normal distributions
      normal_to_lognormal_mean <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        meanlog <- log(normmean ^ 2 / phi)
        return(meanlog)
      }
      
      normal_to_lognormal_sd <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
        return(sdlog)
      }
      
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(theta_tmp$mean[par], theta_tmp$sd[par]), normal_to_lognormal_sd(theta_tmp$mean[par], theta_tmp$sd[par]))
      par_names[par] <- theta_tmp$name[par]
    }
    if(theta_tmp$distribution[par]==""){
      pars[par] <- theta_tmp$mean[par]
      par_names[par] <- theta_tmp$name[par]
    }
    if(theta_tmp$distribution[par]=="normal"){
      pars[par] <- rnorm(1, theta_tmp$mean[par], theta_tmp$sd[par])
      par_names[par] <- theta_tmp$name[par]
    }
  }
  # all_pars <- tibble(i=i, t=t, e=e, ag=ag, name=par_names, value=pars) ## This line doesn't work
  all_pars <- tibble(i=rep(i,nrow(theta_tmp)),t=rep(t,nrow(theta_tmp)), e=rep(e,nrow(theta_tmp)), ag=theta_tmp$antigen_id, name=par_names, value=pars, realized_value=pars) 
  return(all_pars)
}
#' Draw Parameters Fixed Effects With Titer-Dependent Boosting
#'  
#' @description This function draws all parameters listed in theta for successful exposures with fixed effects and titer-dependent boosting.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param ag antigen
#' @param demography Demography information 
#' @param antibody_states An array of true antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_fixed_fx_titer_dep <- function(i, t, e, ag, demography, antibody_states, theta, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  realized <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    pars[par] <- theta_tmp$mean[par]
    realized[par] <- theta_tmp$mean[par]
    par_names[par] <- theta_tmp$name[par]
    
    if(par_names[par] %in% c("boost_short","boost_long")){
      ## Pull out all antigen
      Antigen<-theta_tmp$antigen_id[par]
      titer_threshold <- min(antibody_states[i,t,Antigen], theta_tmp[theta_tmp$name=="titer_ceiling_threshold" & theta_tmp$antigen_id==Antigen, "mean"])
      ## Replace realized titer value for boost parameters
      realized[par] <- pars[par]*(1-theta_tmp[theta_tmp$name=="titer_ceiling_gradient" & theta_tmp$antigen_id==Antigen, "mean"]*titer_threshold)
    }
  }
  all_pars <- tibble(i=rep(i,nrow(theta_tmp)),t=rep(t,nrow(theta_tmp)), e=rep(e,nrow(theta_tmp)), ag=theta_tmp$antigen_id, name=par_names, value=pars, realized_value=realized) 
  return(all_pars)
}

#' Draw Parameters Random Effects With Titer-Dependent Boosting
#'  
#' @description This function draws all parameters listed in theta for successful exposures with random effects and titer-dependent boosting.
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param ag antigen
#' @param demography Demography information 
#' @param antibody_states An array of true antibody titers for all individuals across all time steps and antigens  
#' @param theta Tibble of antibody kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_random_fx_titer_dep <- function(i, t, e, ag, demography, antibody_states, theta, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  realized <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    if(theta_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
      ## Create functions to convert normal distributions to log-normal distributions
      normal_to_lognormal_mean <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        meanlog <- log(normmean ^ 2 / phi)
        return(meanlog)
      }
      
      normal_to_lognormal_sd <- function(normmean, normsd) {
        phi <- sqrt(normsd ^ 2 + normmean ^ 2)
        sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
        return(sdlog)
      }
      
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(theta_tmp$mean[par], theta_tmp$sd[par]), normal_to_lognormal_sd(theta_tmp$mean[par], theta_tmp$sd[par]))
      par_names[par] <- theta_tmp$name[par]
    }
    if(theta_tmp$distribution[par]==""){
      pars[par] <- theta_tmp$mean[par]
      par_names[par] <- theta_tmp$name[par]
    }
    if(theta_tmp$distribution[par]=="normal"){
      pars[par] <- rnorm(1, theta_tmp$mean[par], theta_tmp$sd[par])
      par_names[par] <- theta_tmp$name[par]
    }
    if(par_names[par] %in% c("boost_short","boost_long")){
      ## Pull out all antigen
      Antigen<-theta_tmp$antigen_id[par]
      t1<-t-1
      titer_threshold <- min(antibody_states[i,t1,Antigen], theta_tmp[theta_tmp$name=="titer_ceiling_threshold" & theta_tmp$antigen_id==Antigen, "mean"])
      realized[par] <- pars[par]*(1-theta_tmp[theta_tmp$name=="titer_ceiling_gradient" & theta_tmp$antigen_id==Antigen, "mean"]*titer_threshold)
    }
    if(!(par_names[par] %in% c("boost_short","boost_long"))){
      ## Non boost Realized parameters don't get affected by titer ceiling
      realized[par] <- pars[par]
    }
  }
  all_pars <- tibble(i=rep(i,nrow(theta_tmp)),t=rep(t,nrow(theta_tmp)), e=rep(e,nrow(theta_tmp)), ag=theta_tmp$antigen_id, name=par_names, value=pars, realized_value=realized) 
  return(all_pars)
}

