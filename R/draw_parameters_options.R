#' Draw parameters with fixed effects
#'  
#' @description Draws parameters directly from the distributions specified by `model_pars` for the antibody model with fixed effects. This function ensures that all individuals have the same parameters.
#'
#' @param i individual
#' @param t time
#' @param x exposure
#' @param b biomarker
#' @param demography demography information 
#' @param biomarker_states an array of true biomarker quantities for all individuals across all time steps and biomarkers  
#' @param model_pars tibble of biomarker (antibody) kinetics parameters with variables: 1) exposure_id: numeric exposure ID; 2) biomarker_id: numeric biomarker ID; 3) name: the character name of the parameter; 4) mean: numeric mean of this parameter distribution; 5) sd: the numeric standard deviation of the parameter distribution; 6) distribution: character description of the parameter distribution type (e.g., log-normal, normal)
#' @param ... Additional arguments
#'
#' @return A tibble with the simulated parameters for this exposure event
#' @importFrom dplyr tibble
#' @export
#' @family draw_parameters
#' @examples
#' draw_parameters_fixed_fx(1,1,1,1,example_demography, example_biomarker_states, example_model_pars_numeric)
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
#' Draw parameters with random effects
#'  
#' @description Draws parameters directly from `model_pars` for the antibody model with random effects. Parameters are drawn randomly from a distribution with mean and standard deviation specified within `model_pars`.
#'
#' @inheritParams draw_parameters_fixed_fx
#'
#' @return A tibble with the simulated parameters for this exposure event
#' @importFrom dplyr tibble
#' @importFrom stats rlnorm
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' draw_parameters_random_fx(1,1,1,1,example_demography, example_biomarker_states, example_model_pars_numeric)
draw_parameters_random_fx<- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),]
    pars <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
    if(model_pars_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
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

#' Draw parameters fixed effects with biomarker-quantity-dependent boosting
#'  
#' @description Adds biomarker quantity ceiling effects to the `draw_parameters_fixed_fx` function. Here, an individual’s realized biomarker boost is dependent on their biomarker quantity at the time of the exposure event.
#' 
#' @inheritParams draw_parameters_fixed_fx
#'
#' @return A tibble with the simulated parameters for this exposure event
#' @importFrom dplyr tibble
#' @export
#'
#' @examples
#' model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
#' draw_parameters_fixed_fx_biomarker_dep(2,100,1,1,example_demography, example_biomarker_states_wide, model_pars)
draw_parameters_fixed_fx_biomarker_dep <- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- as.data.frame(model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),])
    
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

#' Draw parameters random effects with biomarker-quantity-dependent boosting
#'  
#' @description Adds biomarker quantity ceiling effects to the previous draw_parameters_random_fx function. Here an individual’s realized biomarker boost is dependent on their biomarker quantity at the time of the exposure event.
#'
#' @inheritParams draw_parameters_fixed_fx
#'
#' @return A tibble with the simulated parameters for this exposure event
#' @importFrom dplyr tibble
#' @importFrom stats rlnorm
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
#' draw_parameters_random_fx_biomarker_dep(2,100,1,1,example_demography, example_biomarker_states_wide, model_pars)
draw_parameters_random_fx_biomarker_dep <- function(i, t, x, b, demography, biomarker_states, model_pars, ...){
  ## Filter for only exposure stimulated 
    model_pars_tmp <- as.data.frame(model_pars[model_pars$exposure_id == x & !is.na(model_pars$exposure_id),])
    pars <- numeric(nrow(model_pars_tmp))
  realized <- numeric(nrow(model_pars_tmp))
  par_names <- character(nrow(model_pars_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(model_pars_tmp)){
      if(is.na(model_pars_tmp$distribution[par]) | model_pars_tmp$distribution[par]==""){
          pars[par] <- model_pars_tmp$mean[par]
          par_names[par] <- model_pars_tmp$name[par]
      } else if(model_pars_tmp$distribution[par] == "log-normal"){ #Convert the normal distributions to log-normal distributions 
      
      ## Create functions to convert normal distributions to log-normal distributions
      pars[par] <- rlnorm(1, normal_to_lognormal_mean(model_pars_tmp$mean[par], model_pars_tmp$sd[par]), normal_to_lognormal_sd(model_pars_tmp$mean[par], model_pars_tmp$sd[par]))
      par_names[par] <- model_pars_tmp$name[par]
    } else if(model_pars_tmp$distribution[par]=="normal"){
      pars[par] <- rnorm(1, model_pars_tmp$mean[par], model_pars_tmp$sd[par])
      par_names[par] <- model_pars_tmp$name[par]
    } else {
        pars[par] <- model_pars_tmp$mean[par]
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

