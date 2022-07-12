#' Simulate Titre Boost (NOT NEEDED ANYMORE)
#'
#' @param i Individual 
#' @param t1 Current time step
#' @param e Exposure_ID
#' @param ag Antigen
#' @param theta Antibody kinetics parameters 
#' @param exposure_histories 
#' @param kinetics_parameters Individual antibody kinetics parameters 
#' @param antigen_map 
#' @param .... 
#'
#' @return
#' @export
#'
#' @examples
simulate_boost <- function(i,t1, e, ag, theta, exposure_histories,kinetics_parameters,antigen_map, ....){
  boost <- draw_parameters(i, t1, e, ag, demography, theta, antibody_state, ...)
  titre_threshold <- min(antibody_states[i,t1,ag], theta[theta$name=="titre_ceiling_threshold" & theta$e==e & theta$ag==ag, "mean"])
  boost <- boost*(1-theta[theta$name=="titre_ceiling_gradient" & theta$e==e & theta$ag==ag, "mean"]*titre_threshold)
  boost
}






