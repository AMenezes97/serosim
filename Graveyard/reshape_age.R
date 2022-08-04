#' Reshape Age Matrix
#'
#' @param age A matrix containing an individual's age at each time step
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return An age data set is returned
#' @export
#'
#' @examples
reshape_age <-function(age, N_pathogens){
  for(pathogen in 1:N_pathogens){
    age[[pathogen]] <- reshape2::melt(age[[pathogen]])
    colnames(age[[pathogen]]) <- c("Individual","Time","Age")
    age[[pathogen]]$Pathogen <- pathogen
    age[[pathogen]]$Age <- as.factor(age[[pathogen]]$Age) #check to make sure this works without adding 1 to all ages
  }

  age <- do.call("bind_rows",age)
  age$Pathogen <- as.factor(age$Pathogen)
  return(age)
}
