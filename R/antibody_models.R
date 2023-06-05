#' Monophasic antibody boosting-waning model
#' 
#' @description This monophasic antibody boosting-waning model model assumes that for each exposure there is a boost and boost waning parameter
#'
#' @param i individual
#' @param t1 time
#' @param b biomarker
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param biomarker_states An array of biomarker states (biomarker quantities) across all individuals, time steps and biomarker IDs
#' @param kinetics_parameters A tibble of parameters needed for the antibody kinetics model for all biomarkers 
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param ... Additional arguments
#'
#' @return A biomarker quantity is returned 
#' @export
#'
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx(1,1,1,1,NULL, NULL, example_model_pars_numeric)
#' antibody_model_monophasic(1,1,1,example_exposure_histories_wide, example_biomarker_states_wide, 
#' tmp_pars, example_biomarker_map_numeric)
antibody_model_monophasic <-  function(i, t1, b, exposure_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){
  ## Find which successful exposures correspond to this biomarker 
  exposure_id_tmp<-biomarker_map$exposure_id[biomarker_map$biomarker_id==b]
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting biomarker quantity to 0
  biomarker_quantity<-0
  
  ## Calculate current biomarker quantity if there has been an exposure 
  if(sum(exp_history,na.rm = TRUE)==0){
    return(0)
  }
  if(sum(exp_history,na.rm = TRUE)>0){
    ## Extract all kinetics_parameters for biomarker 
    b_tmp<-b
    
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
    tmp_kinetics_parameters<-tmp_kinetics_parameters[tmp_kinetics_parameters$b==b_tmp,] 
    
    tmp_boost <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",] 
    tmp_wane <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane",] 
    
    for(j in seq_along(tmp_boost$realized_value)){
      biomarker_quantity<- biomarker_quantity + tmp_boost$realized_value[j]*max(0,1-tmp_wane$realized_value[j]*(t1-tmp_wane$t[j]))
    }
    biomarker_quantity
  }
}

#' Biphasic antibody boosting-waning model
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @inheritParams antibody_model_monophasic
#' 
#' 
#' @return A biomarker quantity is returned 
#' @importFrom data.table data.table
#' @export
#'
#' @examples
#' model_pars <- reformat_biomarker_map(example_model_pars_biphasic)
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_fixed_fx_biomarker_dep(1,1,1,1,NULL, NULL, model_pars)
#' antibody_model_biphasic(1,1,1,example_exposure_histories_wide, example_biomarker_states_wide, 
#' tmp_pars, example_biomarker_map_numeric)
antibody_model_biphasic <-  function(i, t1, b, exposure_histories, biomarker_states, kinetics_parameters, biomarker_map, ...){

  ## Find which successful exposures correspond to this biomarker 
  exposure_id_tmp<-biomarker_map$exposure_id[biomarker_map$biomarker_id==b]
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting biomarker quantity to 0
  biomarker_quantity<-0
  
  ## Calculate current biomarker quantity if there has been an exposure 
  if(sum(exp_history,na.rm = TRUE)==0){
    return(0)
  }
  if(sum(exp_history,na.rm = TRUE)>0){
    ## Extract all kinetics_parameters for biomarker 
    b_tmp<-b
    
    tmp_kinetics_parameters <- data.table(kinetics_parameters[[i]])
    tmp_kinetics_parameters<-tmp_kinetics_parameters[tmp_kinetics_parameters$b==b_tmp,] 
    
    tmp_boost_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",] 
    tmp_boost_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_short",] 
    
    tmp_wane_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_long",] 
    tmp_wane_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_short",] 
    
    for(j in seq_along(tmp_boost_long$realized_value)){
      biomarker_quantity<- biomarker_quantity + tmp_boost_long$realized_value[j]*max(0,1-tmp_wane_long$realized_value[j]*(t1-tmp_wane_long$t[j])) + tmp_boost_short$realized_value[j]*max(0,1-tmp_wane_short$realized_value[j]*(t1-tmp_wane_short$t[j]))
    }
    biomarker_quantity
    
  }
}




#' Power law antibody boosting and waning
#' 
#' @description Power law antibody boosting and waning model
#'
#' @inheritParams antibody_model_monophasic
#' @param t time 
#'
#' @return A biomarker_quantity is returned 
#' @export
#'
#' @examples
#' tmp_pars <- list()
#' tmp_pars[[1]] <- draw_parameters_random_fx(1,1,1,1,NULL,NULL,example_model_pars_typhoid)
#' tmp_exposure_history <- array(0,dim=c(1,11,2))
#' tmp_exposure_history[1,1,1] <- 1
#' antibody_model_typhoid(1,10, 1, tmp_exposure_history, NULL, tmp_pars,example_biomarker_map_numeric)
antibody_model_typhoid <- function(i, t, b, exposure_histories=NULL, biomarker_states=NULL, kinetics_parameters, biomarker_map=NULL,...){
    tmp_pars <- kinetics_parameters[[i]]
    
    ## Find which successful exposures correspond to this biomarker 
    exposure_id_tmp<-biomarker_map$exposure_id[biomarker_map$biomarker_id==b]
    
    ## Find all exposures up until current time for this individual and exposure type
    exp_history <- exposure_histories[i,1:t,exposure_id_tmp]
    
    if(sum(exp_history,na.rm = TRUE)==0){
        return(0)
    }
    
    if(nrow(tmp_pars) > 1){
        tmp_pars <- as.data.frame(tmp_pars)
        
        titer <- tmp_pars[tmp_pars$b == b & tmp_pars$name == "y0","value"][1]
        ## There will be a distinct set of parameters for each exposure in exposure history
        ## Get exposure parameters relevant to this biomarker
        tmp_pars <- tmp_pars[tmp_pars$b == b & tmp_pars$t <= t,]
        
        ## Assume that tmp_pars is in the correct time order
        ##########
        ## Unlike other models, the time order of this one matters.
        ## We just use the most recent boosting event to calculate the current titer. But this
        ## means that we also need to know the titer at the time of the boost
        y1s <- tmp_pars[tmp_pars$name == "y1","value"]
        alphas <- tmp_pars[tmp_pars$name == "alpha","value"]
        rs <- tmp_pars[tmp_pars$name == "r","value"]
        t1s <- tmp_pars[tmp_pars$name == "t1","value"]
        t_inf <- tmp_pars[tmp_pars$name == "y1","t"]
        typhoid <- function(t, y0, y1, beta, r, t1){
            mu <- (1/t1) * log(y1/y0)
            #y <- numeric(length(t))
            #y[t<=t1] <- y0*exp(mu*t[t <= t1])
            alpha <- 1/(r-1)
            if(t <= t1){
                y <- y0*exp(mu*t)
            } else {
                tau <- t - t1
                y <- y1*(1+beta*tau)^(-alpha)
            }
            y
        }
        y0 <- as.numeric(titer)
        if(length(t_inf)>1){
            for(x in 2:length(t_inf)){
                y0 <- typhoid(t_inf[x]-t_inf[x-1], y0, y1s[x-1], alphas[x-1], rs[x-1], t1s[x-1])
            }
        }
        y <- typhoid(t-t_inf[length(t_inf)], y0, y1s[length(t1s)], alphas[length(alphas)],rs[length(rs)],t1s[length(t1s)])
    } else {
        y <- titer
    }
    as.numeric(y)
    
}
