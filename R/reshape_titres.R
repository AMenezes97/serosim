#' Reshape Titres Matrix
#'
#' @param titres A matrix containing the antibody titres for each individual at each time step
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return A titre data set is returned
#' @export
#'
#' @examples
reshape_titres <-function(titres, N_pathogens){
  for(pathogen in 1:N_pathogens){
    titres[[pathogen]] <- reshape2::melt(titres[[pathogen]])
    colnames(titres[[pathogen]]) <- c("Individual","Time","Titre")
    titres[[pathogen]]$Pathogen <- pathogen
  }

  titres <- do.call("bind_rows",titres)
  titres$Pathogen <- as.factor(titres$Pathogen)
  return(titres)
}
