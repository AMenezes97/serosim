#' Reshape Infection Histories Across Time For All Individuals And Pathogens
#'
#' @param infection_histories A matrix containing the infection histories for each individual at each time step
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return An infection history data set is returned
#' @export
#'
#' @examples
reshape_infection_histories <- function(infection_histories, N_pathogens){
  for(pathogen in 1:N_pathogens){
    infection_histories[[pathogen]] <- reshape2::melt(infection_histories[[pathogen]])
    colnames(infection_histories[[pathogen]]) <- c("Individual","Time","Infected?")
    infection_histories[[pathogen]]$Pathogen <- pathogen
    infection_histories[[pathogen]]$`Infected?` <- as.factor(infection_histories[[pathogen]]$`Infected?`)
  }

  infection_histories <- do.call("bind_rows",infection_histories)
  infection_histories$Pathogen <- as.factor(infection_histories$Pathogen)
  return(infection_histories)
}
