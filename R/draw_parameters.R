#' Draw Antibody Parameters
#'
#' @param i Individual 
#' @param t Time
#' @param e Exposure_ID
#' @param ag Antigen
#' @param demography Population demography dat set 
#' @param theta Antibody kinetics parameters
#' @param antibody_state 
#' @param ... 
#'
#' @return 
#' @export
#'
#' @examples 
draw_parameters <- function(i, t, e, ag, demography, theta, antibody_state, ...){
  ## Filter for only exposure stimulated 
  theta_tmp <- theta %>% filter(exposure_id == e)
  pars <- numeric(nrow(theta_tmp))
  par_names <- character(nrow(theta_tmp))
  ## For each parameter; randomly sample from the distribution given the mean and sd 
  for(par in 1:nrow(theta_tmp)){
    if(theta_tmp$distribution == "log-normal"){
      pars[par] <- rlnorm(1, theta_tmp$mean[par], theta_tmp$sd[par])
      par_names[par] <- theta_tmp$name[par]
    }
    else{
      pars[par] <- rnorm(1, theta_tmp$mean[par], theta_tmp$sd[par])
      par_names[par] <- theta_tmp$name[par]
    }
    #Add titre-dependent boosting  
      if(par_names[par] %in% c("boost_short","boost_long")){
        titre_threshold <- min(antibody_states[i,t1,ag], theta_tmp[theta_tmp$name=="titre_ceiling_threshold" & theta_tmp$ag==ag, "mean"])
        pars[par] <- pars[par]*(1-theta_tmp[theta_tmp$name=="titre_ceiling_gradient" & theta_tmp$ag==ag, "mean"]*titre_threshold)
    }
  }
  all_pars <- tibble(i=i, t=t, e=e, ag=ag, name=par_names, value=pars)
  return(all_pars)
}
