#' Draw Parameters Fixed Effects
#'  
#' @description This function draws parameters directly from model_pars for the antibody model with fixed effects. This function ensures that all individuals have the same parameters.
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param b biomarker
#' @param demography Demography information 
#' @param biomarker_states An array of true biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of biomarker (antibody) kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_fixed_fx <- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
  model_pars_tmp <- model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),]
  pars <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
    pars[par] <- model_pars_tmp$mean[par]
    par_names[par] <- model_pars_tmp$name[par]
  }
  all_pars <- tibble(i=rep(i,nrow(model_pars_tmp)),t=rep(t,nrow(model_pars_tmp)), x=rep(x,nrow(model_pars_tmp)), b=model_pars_tmp$biomarker_id, name=par_names, value=pars, realized_value=pars) 
  return(all_pars)
}
#' Draw Parameters Random Effects
#'  
#' @description This function draws parameters directly from model_pars for the antibody model with random effects. Parameters are drawn randomly from a distribution with mean and standard deviation specified within model_pars.
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param b biomarker
#' @param demography Demography information 
#' @param biomarker_states An array of true biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of biomarker (antibody) kinetics parameters  
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_random_fx<- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),]
    pars <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
    if(model_pars_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
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
      
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(model_pars_tmp$mean[par], model_pars_tmp$sd[par]), normal_to_lognormal_sd(model_pars_tmp$mean[par], model_pars_tmp$sd[par]))
      par_names[par] <- model_pars_tmp$name[par]
    }
    if(model_pars_tmp$distribution[par]==""){
      pars[par] <- model_pars_tmp$mean[par]
      par_names[par] <- model_pars_tmp$name[par]
    }
    if(model_pars_tmp$distribution[par]=="normal"){
      pars[par] <- rnorm(1, model_pars_tmp$mean[par], model_pars_tmp$sd[par])
      par_names[par] <- model_pars_tmp$name[par]
    }
  }
  all_pars <- tibble(i=rep(i,nrow(model_pars_tmp)),t=rep(t,nrow(model_pars_tmp)), x=rep(x,nrow(model_pars_tmp)), b=model_pars_tmp$biomarker_id, name=par_names, value=pars, realized_value=pars) 
  return(all_pars)
}

#' Draw Parameters Fixed Effects With Biomarker Quantity Dependent Boosting
#'  
#' @description This function adds biomarker quantity ceiling effects to the previous draw_parameters_fixed_fx function. Here an individual’s realized biomarker boost is dependent on their biomarker quantity at the time of the exposure event.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param b biomarker
#' @param demography Demography information 
#' @param biomarker_states An array of true biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of biomarker (antibody) kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_fixed_fx_biomarker_dep <- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),]
    
  pars <- numeric(nrow(model_pars_tmp))
  realized <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
    pars[par] <- model_pars_tmp$mean[par]
    realized[par] <- model_pars_tmp$mean[par]
    par_names[par] <- model_pars_tmp$name[par]
    
    if(par_names[par] %in% c("boost_short","boost_long","boost")){
      ## Pull out all biomarker
      biomarker<-model_pars_tmp$biomarker_id[par]
      biomarker_threshold <- min(biomarker_states[i,t,biomarker], model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_threshold" & model_pars_tmp$biomarker_id==biomarker, "mean"])
      ## Replace realized biomarker quantity for boost parameters
      realized[par] <- pars[par]*(1-model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_gradient" & model_pars_tmp$biomarker_id==biomarker, "mean"]*biomarker_threshold)
    }
  }
  all_pars <- tibble(i=rep(i,nrow(model_pars_tmp)),t=rep(t,nrow(model_pars_tmp)), x=rep(x,nrow(model_pars_tmp)), b=model_pars_tmp$biomarker_id, name=par_names, value=pars, realized_value=realized) 
  return(all_pars)
}

#' Draw Parameters Random Effects With Biomarker Quantity Dependent Boosting
#'  
#' @description his function adds biomarker quantity ceiling effects to the previous draw_parameters_random_fx function. Here an individual’s realized biomarker boost is dependent on their biomarker quantity  at the time of the exposure event.
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param b biomarker
#' @param demography Demography information 
#' @param biomarker_states An array of true biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars Tibble of biomarker (antibody) kinetics parameters 
#' @param ... 
#'
#' @return A tibble with the parameters drawn is returned
#' @export
#'
#' @examples
draw_parameters_random_fx_biomarker_dep <- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),]
      pars <- numeric(nrow(model_pars_tmp))
  realized <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
    if(model_pars_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
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
      
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(model_pars_tmp$mean[par], model_pars_tmp$sd[par]), normal_to_lognormal_sd(model_pars_tmp$mean[par], model_pars_tmp$sd[par]))
      par_names[par] <- model_pars_tmp$name[par]
    }
    if(model_pars_tmp$distribution[par]==""){
      pars[par] <- model_pars_tmp$mean[par]
      par_names[par] <- model_pars_tmp$name[par]
    }
    if(model_pars_tmp$distribution[par]=="normal"){
      pars[par] <- rnorm(1, model_pars_tmp$mean[par], model_pars_tmp$sd[par])
      par_names[par] <- model_pars_tmp$name[par]
    }
    if(par_names[par] %in% c("boost_short","boost_long","boost")){
      ## Pull out all biomarker
      biomarker<-model_pars_tmp$biomarker_id[par]
      t1<-t-1
      biomarker_threshold <- min(biomarker_states[i,t1,biomarker], model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_threshold" & model_pars_tmp$biomarker_id==biomarker, "mean"])
      realized[par] <- pars[par]*(1-model_pars_tmp[model_pars_tmp$name=="biomarker_ceiling_gradient" & model_pars_tmp$biomarker_id==biomarker, "mean"]*biomarker_threshold)
    }
    if(!(par_names[par] %in% c("boost_short","boost_long","boost"))){
      ## Non-boost realized parameters don't get affected by biomarker quantity ceiling
      realized[par] <- pars[par]
    }
  }
  all_pars <- tibble(i=rep(i,nrow(model_pars_tmp)),t=rep(t,nrow(model_pars_tmp)), x=rep(x,nrow(model_pars_tmp)), b=model_pars_tmp$biomarker_id, name=par_names, value=pars, realized_value=realized) 
  return(all_pars)
}

