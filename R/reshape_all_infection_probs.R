#' Reshape All Infection Probabilities Matrix
#'
#' @param all_infection_probs A matrix containing the probability of infection for each individual at each time step
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return An infection probabilities data set is returned
#' @export
#'
#' @examples
reshape_all_infection_probs <- function(all_infection_probs, N_pathogens){
  for(pathogen in 1:N_pathogens){
    all_infection_probs[[pathogen]] <- reshape2::melt(all_infection_probs[[pathogen]])
    colnames(all_infection_probs[[pathogen]]) <- c("Individual","Time","Infection probability")
    all_infection_probs[[pathogen]]$Pathogen <- pathogen
  }
  all_infection_probs <- do.call("bind_rows",all_infection_probs)
  all_infection_probs$Pathogen <- as.factor(all_infection_probs$Pathogen)
return(all_infection_probs)
}
