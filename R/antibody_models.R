#' Monophasic antibody boosting-waning model
#' 
#' @description This monophasic antibody boosting-waning model assumes that for each exposure there is a boost and waning parameter
#'
#' @param i individual
#' @param t1 time
#' @param b biomarker
#' @param immune_histories An array of immune histories across all individuals, time steps and exposure IDs
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs
#' @param kinetics_parameters A tibble of parameters needed for the antibody kinetics model for all biomarkers 
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param ... Additional arguments
#'
#' @return A biomarker quantity is returned 
#' @seealso antibody_model_monophasic_cpp
#' @family antibody_models
#' @export
#'
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, example_model_pars_numeric)
#' antibody_model_monophasic(1,1,1,example_immune_histories_wide, example_biomarker_states_wide, 
#' tmp_pars, example_biomarker_map_numeric)
antibody_model_monophasic <-  function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  biomarker_quantity <- 0
  
  ## Get kinetics parameters for this individual
  if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
  } else {
    return(biomarker_quantity)
  }
  
  ## Get only exposures relevant to this biomarker ID and time
  tmp_kinetics_parameters <- tmp_kinetics_parameters[tmp_kinetics_parameters$b == b & tmp_kinetics_parameters$t <= t1,]
  
  ## Only continue if there are relevant exposures to calculate kinetics for
  if(nrow(tmp_kinetics_parameters) > 0){
    boosts <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",]$realized_value ## Boosts
    wanes <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane",]$realized_value ## Waning rate
    t_infs <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane",]$t ## Time of infections
    ## Sum up contribution of each boost, with waning
    for(j in seq_along(t_infs)){
      biomarker_quantity<- biomarker_quantity + boosts[j]*max(0,1-wanes[j]*(t1-t_infs[j]))
    }
  }
  biomarker_quantity
}

#' Biphasic antibody boosting-waning model
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term waning, short-term boost, and short-term waning parameters
#'
#' @inheritParams antibody_model_monophasic
#' 
#' 
#' @return A biomarker quantity is returned 
#' @importFrom data.table data.table
#' @family antibody_models
#' @export
#'
#' @examples
#' model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx_biomarker_dep(1,1,1,NULL, NULL,model_pars)
#' antibody_model_biphasic(1,1,1,example_immune_histories_wide, example_biomarker_states_wide,
#' tmp_pars, model_pars)
antibody_model_biphasic <-  function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  biomarker_quantity <- 0
  
  ## Get kinetics parameters for this individual
  if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
  } else {
    return(biomarker_quantity)
  }
  
  ## Get only exposures relevant to this biomarker ID and time
  tmp_kinetics_parameters <- tmp_kinetics_parameters[tmp_kinetics_parameters$b == b & tmp_kinetics_parameters$t <= t1,]
  
  ## Only continue if there are relevant exposures to calculate kinetics for
  if(nrow(tmp_kinetics_parameters) > 0){
    boosts_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",]$realized_value ## Boosts
    boosts_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_short",]$realized_value ## Boosts
    
    wanes_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_long",]$realized_value ## Waning rate
    wanes_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_short",]$realized_value ## Waning rate
    
    t_infs <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",]$t ## Time of infections
    ## Sum up contribution of each boost, with waning
    for(j in seq_along(t_infs)){
      biomarker_quantity<- biomarker_quantity + boosts_long[j]*max(0,1-wanes_long[j]*(t1-t_infs[j])) + boosts_short[j]*max(0,1-wanes_short[j]*(t1-t_infs[j]))
    }
  }
  biomarker_quantity
}


typhoid <- function(t, y0, y1, beta, r, t1){
  mu <- (1/t1) * log(y1/y0)
  alpha <- 1/(r-1)
  if(t <= t1){
    y <- y0*exp(mu*t)
  } else {
    tau <- t - t1
    y <- y1*(1+beta*tau)^(-alpha)
  }
  y
}

#' Power law antibody boosting and waning
#' 
#' @description Power law antibody boosting and waning model
#'
#' @inheritParams antibody_model_monophasic
#'
#' @return A biomarker_quantity is returned 
#' @export
#' @family antibody_models
#'
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_random_fx(1,1,1,NULL,NULL,example_model_pars_typhoid)
#' tmp_immune_history <- array(0,dim=c(1,11,2))
#' tmp_immune_history[1,1,1] <- 1
#' antibody_model_typhoid(1,10, 1, tmp_immune_history, NULL, tmp_pars,example_biomarker_map_numeric)
antibody_model_typhoid <- function(i, t1, b, immune_histories=NULL, biomarker_states=NULL, kinetics_parameters, biomarker_map=NULL,...){
  
  biomarker_quantity <- 0
  ## Get only exposures relevant to this biomarker ID and time
  if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
  } else {
    return(biomarker_quantity)
  }
  tmp_pars <- tmp_kinetics_parameters[tmp_kinetics_parameters$b == b & tmp_kinetics_parameters$t <= t1,]
  
  ## Only continue if there are relevant exposures to calculate kinetics for
  if(nrow(tmp_pars) > 0){
        titer <- tmp_pars[tmp_pars$b == b & tmp_pars$name == "y0",]$value[1]
        ## Assume that tmp_pars is in the correct time order
        ##########
        ## Unlike other models, the time order of this one matters.
        ## We just use the most recent boosting event to calculate the current titer. But this
        ## means that we also need to know the titer at the time of the boost
        y1s <- tmp_pars[tmp_pars$name == "y1",]$value
        alphas <- tmp_pars[tmp_pars$name == "alpha",]$value
        rs <- tmp_pars[tmp_pars$name == "r",]$value
        t1s <- tmp_pars[tmp_pars$name == "t1",]$value
        t_inf <- tmp_pars[tmp_pars$name == "y1",]$t
        y0 <- as.numeric(titer)
        if(length(t_inf)>1){
            for(x in 2:length(t_inf)){
                y0 <- typhoid(t_inf[x]-t_inf[x-1], y0, y1s[x-1], alphas[x-1], rs[x-1], t1s[x-1])
            }
        }
        biomarker_quantity <- typhoid(t1-t_inf[length(t_inf)], y0, y1s[length(t1s)], alphas[length(alphas)],rs[length(rs)],t1s[length(t1s)])
  }
  biomarker_quantity
}

#' Monophasic antibody boosting-waning model with cross-reactive strains
#' 
#' @description Monophasic antibody boosting-waning model with cross-reactivity between strains. This monophasic antibody boosting-waning model assumes that for each exposure there is a boost and waning parameter describing antibody kinetics against the infecting strain (i.e., for exposure_id==biomarker_id). The model loops through each exposure type and reduces the amount of boosting as a function of cross-reactivity, which is determined by a proportion given in the `biomarker_map` data frame as the `value` variable.
#'
#' @inheritParams antibody_model_monophasic
#' 
#' 
#' @return A biomarker quantity is returned 
#' @importFrom data.table data.table
#' @export
#' @family antibody_models
#'
#' @examples
#' tmp_pars <- list()
#' ## Set up simple model_pars table for this antibody model
#' library(dplyr)
#' model_pars_tmp <- example_model_pars_numeric %>% mutate(biomarker_id = exposure_id)
#' ## Simulate one infection with exposure ID 1 at t=1
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, model_pars_tmp)
#'  
#' ## Set up a simple biomarker map for cross-reactivity
#' biomarker_map = tidyr::expand_grid(exposure_id=1:2, biomarker_id=1:2) %>% 
#' mutate(value = if_else(exposure_id==biomarker_id, 1, 0.5))
#' antibody_model_monophasic_cross_reactivity(1,1,1,
#' example_immune_histories_wide, example_biomarker_states_wide, 
#' tmp_pars, biomarker_map)
antibody_model_monophasic_cross_reactivity <-  function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  biomarker_quantity <- 0
  
  ## Get kinetics parameters for this individual
  if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
  } else {
    return(biomarker_quantity)
  }

  ## Only use exposure IDs relevant to this biomarker. In this model, this should give entries where b==x
  use_exposure_ids <- unique(biomarker_map[biomarker_map$biomarker_id == b, ]$exposure_id)
    
  ## Get only exposures relevant to this biomarker ID and time
  tmp_kinetics_parameters <- tmp_kinetics_parameters[tmp_kinetics_parameters$x %in% use_exposure_ids & tmp_kinetics_parameters$t <= t1,]
  ## Only continue if there are relevant exposures to calculate kinetics for
  if(nrow(tmp_kinetics_parameters) > 0){
    boosts <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",]$realized_value ## Boosts
    wanes <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane",]$realized_value ## Waning rate
    t_infs <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",]$t ## Time of infections
    exposure_ids <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",]$x ## Exposure ID of infections
    
    ## Sum up contribution of each boost, with waning
    for(j in seq_along(t_infs)){
      ## In this model, cross-reactivity is just the homologous boost/wane amount multiplied by a proportion
      cross_reactivity <- biomarker_map[biomarker_map$exposure_id == exposure_ids[j] & 
                                          biomarker_map$biomarker_id == b,]$value
      biomarker_quantity<- biomarker_quantity + max(0,boosts[j]*cross_reactivity-wanes[j]*(t1-t_infs[j]))
    }
  }
  biomarker_quantity
}

#' Biphasic antibody boosting-waning model with cross-reactive strains
#' 
#' @description Biphasic antibody boosting-waning model with cross-reactivity between strains. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters describing antibody kinetics against the infecting strain (i.e., for exposure_id==biomarker_id). The model loops through each exposure type and reduces the amount of boosting as a function of cross-reactivity, which is determined by a proportion given in the `biomarker_map` data frame as the `value` variable.
#'
#' @inheritParams antibody_model_monophasic
#' 
#' 
#' @return A biomarker quantity is returned 
#' @importFrom data.table data.table
#' @export
#' @family antibody_models
#'
#' @examples
#' tmp_pars <- list()
#' ## Set up simple model_pars table for this antibody model
#' library(dplyr)
#' model_pars_tmp <- example_model_pars_biphasic %>% reformat_biomarker_map() %>% 
#' mutate(biomarker_id = exposure_id)
#' ## Simulate one infection with exposure ID 1 at t=1
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, model_pars_tmp)
#'  
#' ## Set up a simple biomarker map for cross-reactivity
#' biomarker_map = tidyr::expand_grid(exposure_id=1:2, biomarker_id=1:2) %>% 
#' mutate(value = if_else(exposure_id==biomarker_id, 1, 0.5))
#' antibody_model_biphasic_cross_reactivity(1,1,1,example_immune_histories_wide, 
#' example_biomarker_states_wide, tmp_pars, biomarker_map)
antibody_model_biphasic_cross_reactivity <-  function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  biomarker_quantity <- 0
  
  ## Get kinetics parameters for this individual
  if(!is.null(kinetics_parameters[[i]])){
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
  } else {
    return(biomarker_quantity)
  }
  
  ## Only use exposure IDs relevant to this biomarker. In this model, this should give entries where b==x
  use_exposure_ids <- unique(biomarker_map[biomarker_map$biomarker_id == b, ]$exposure_id)
  
  ## Get only exposures relevant to this biomarker ID and time
  tmp_kinetics_parameters <- tmp_kinetics_parameters[tmp_kinetics_parameters$x %in% use_exposure_ids & tmp_kinetics_parameters$t <= t1,]
  ## Only continue if there are relevant exposures to calculate kinetics for
  if(nrow(tmp_kinetics_parameters) > 0){
    boosts_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",]$realized_value ## Boosts
    boosts_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_short",]$realized_value ## Boosts
    
    wanes_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_long",]$realized_value ## Waning rate
    wanes_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_short",]$realized_value ## Waning rate
    
    t_infs <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",]$t ## Time of infections
    
    exposure_ids <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",]$x ## Exposure ID of infections
    
    ## Sum up contribution of each boost, with waning
    for(j in seq_along(t_infs)){
      ## In this model, cross-reactivity is just the homologous boost/wane amount multiplied by a proportion
      cross_reactivity <- biomarker_map[biomarker_map$exposure_id == exposure_ids[j] & 
                                          biomarker_map$biomarker_id == b,]$value
      biomarker_quantity<- biomarker_quantity + max(0,boosts_long[j]*cross_reactivity-wanes_long[j]*(t1-t_infs[j])) + max(0,boosts_short[j]*cross_reactivity-wanes_short[j]*(t1-t_infs[j]))
    }
  }
  biomarker_quantity
}

#' Monophasic antibody boosting-waning model Rcpp implementation
#' 
#' @description Identical to \code{\link{antibody_model_monophasic}}, but implemented in Cpp
#'
#' @inheritParams antibody_model_monophasic
#'
#' @return Biomarker quantity at the specified time point
#' @export
#' @family antibody_models
#'
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, example_model_pars_numeric)
#' antibody_model_monophasic_cpp(1,1,1,example_immune_histories_wide, example_biomarker_states_wide,tmp_pars, example_biomarker_map_numeric)
antibody_model_monophasic_cpp <- function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  antibody_model_monophasic_cpp_internal(i,t1,b,immune_histories,biomarker_states,kinetics_parameters,biomarker_map)
}


#' Biphasic antibody boosting-waning model Rcpp implementation
#' 
#' @description Identical to \code{\link{antibody_model_biphasic}}, but implemented in Cpp
#'
#' @inheritParams antibody_model_monophasic
#'
#' @return Biomarker quantity at the specified time point
#' @family antibody_models
#' @export
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,NULL, NULL, example_model_pars_numeric)
#' antibody_model_biphasic_cpp(1,1,1,example_immune_histories_wide, example_biomarker_states_wide,tmp_pars, example_biomarker_map_numeric)
antibody_model_biphasic_cpp <- function(i, t1, b, immune_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  antibody_model_biphasic_cpp_internal(i,t1,b,immune_histories,biomarker_states,kinetics_parameters,biomarker_map)
}

