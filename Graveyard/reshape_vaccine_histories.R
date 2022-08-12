#' Reshape Vaccine Histories Across Time For All Individuals And Pathogens
#'
#' @param vaccine_histories A matrix containing the vaccination histories for each individual at each time step
#'
#' @return A vaccine history data set is returned
#' @export
#'
#' @examples
reshape_vaccine_histories <- function(vaccine_histories){
  vaccine_histories_reshaped <- reshape2::melt(vaccine_histories) %>% dplyr::mutate(value=as.factor(value))
  colnames(vaccine_histories_reshaped) <- c("Individual","Time","Vaccinated?")
  return(vaccine_histories_reshaped)
}

