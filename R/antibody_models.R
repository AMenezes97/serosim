#' Monophasic antibody boosting-waning model
#' 
#' @description This monophasic antibody boosting-waning model model assumes that for each exposure there is a boost and boost waning parameter
#'
#' @param i Individual
#' @param t1 time
#' @param b biomarker
#' @param exposure_histories An array of exposure histories across all individuals, time steps and exposure IDs
#' @param antibody_states An array of antibody states across all individuals, time steps and biomarker IDs
#' @param kinetics_parameters A tibble of parameters needed for the antibody kinetics model for all biomarkers 
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param ... 
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_monophasic <-  function(i, t1, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map, ...){
  ## Find which successful exposures correspond to this biomarker 
  exposure_id_tmp<-biomarker_map$exposure_id[biomarker_map$biomarker_id==b]
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting titer to 0
  titer<-0
  
  ## Calculate current titer if there has been an exposure 
  if(sum(exp_history,na.rm = TRUE)==0){
    return(0)
  }
  if(sum(exp_history,na.rm = TRUE)>0){
    ## Extract all kinetics_parameters for biomarker 
    b_tmp<-b
    
    tmp_kinetics_parameters <- kinetics_parameters[[i]]
    tmp_kinetics_parameters<-tmp_kinetics_parameters[tmp_kinetics_parameters$b==b_tmp,] ## Since you are going through time, all parameters will only be from the current or previous times?
    
    # setkey(tmp_kinetics_parameters, cols="i","t","e","b","name","value", "realized_value")
    tmp_boost <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost",] 
    tmp_wane <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane",] 
    
    for(j in seq_along(tmp_boost$realized_value)){
      titer<- titer + tmp_boost$realized_value[j]*max(0,1-tmp_wane$realized_value[j]*(t1-tmp_wane$t[j]))
    }
    titer
  }
}

#' Biphasic antibody boosting-waning model
#' 
#' @description Biphasic antibody boosting-waning model. This model assumes that for each exposure there is a set of long-term boost, long-term boost waning, short-term boost, and short-term boost waning parameters
#'
#' @inheritParams antibody_model_monophasic
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_biphasic <-  function(i, t1, b, exposure_histories, antibody_states, kinetics_parameters, biomarker_map, ...){

  ## Find which successful exposures correspond to this biomarker 
  exposure_id_tmp<-biomarker_map$exposure_id[biomarker_map$biomarker_id==b]
  
  ## Find all exposures up until current time for this individual and exposure type
  exp_history <- exposure_histories[i,1:t1,exposure_id_tmp]
  
  ## Set starting titer to 0
  titer<-0
  
  ## Calculate current titer if there has been an exposure 
  if(sum(exp_history,na.rm = TRUE)==0){
    return(0)
  }
  if(sum(exp_history,na.rm = TRUE)>0){
    ## Extract all kinetics_parameters for biomarker 
    b_tmp<-b
    
    tmp_kinetics_parameters <- data.table(kinetics_parameters[[i]])
    tmp_kinetics_parameters<-tmp_kinetics_parameters[tmp_kinetics_parameters$b==b_tmp,] ## Since you are going through time, all parameters will only be from the current or previous times?
    
    # setkey(tmp_kinetics_parameters, cols="i","t","e","b","name","value", "realized_value")
    tmp_boost_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_long",] 
    tmp_boost_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "boost_short",] 
    
    tmp_wane_long <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_long",] 
    tmp_wane_short <- tmp_kinetics_parameters[tmp_kinetics_parameters$name == "wane_short",] 
    
    
    for(j in seq_along(tmp_boost_long$realized_value)){
      titer<- titer + tmp_boost_long$realized_value[j]*max(0,1-tmp_wane_long$realized_value[j]*(t1-tmp_wane_long$t[j])) + tmp_boost_short$realized_value[j]*max(0,1-tmp_wane_short$realized_value[j]*(t1-tmp_wane_short$t[j]))
    }
    titer
    
  }
}


#' Power law antibody boosting and waning
#' 
#' @description 
#'
#' @inheritParams antibody_model_monophasic
#'
#' @return A titer value is returned 
#' @export
#'
#' @examples
antibody_model_typhoid <- function(i, t, b, exposure_histories=NULL, antibody_states=NULL, kinetics_parameters, biomarker_map=NULL,...){
    tmp_pars <- kinetics_parameters[[i]]
    titer <- tmp_pars[tmp_pars$b == b & tmp_pars$name == "y0","value"] 
    ## There will be a distinct set of parameters for each exposure in exposure history
    ## Get exposure parameters relevant to this biomarker
    tmp_pars <- tmp_pars[tmp_pars$b == b & tmp_pars$t <= t,]
    
    if(nrow(tmp_pars) > 1){
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
            y <- numeric(length(t))
            y[t<=t1] <- y0*exp(mu*t[t <= t1])
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
