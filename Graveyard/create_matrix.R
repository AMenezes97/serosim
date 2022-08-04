#' Create Matrices To Store Simulation Information
#'
#' @param N_pathogens The number of pathogens in the simulation
#' @param N The number of individuals in the simulation
#' @param times The number of time steps in the simulation
#' @param fill Item to fill in matrix with; defaults to NA
#'
#' @return An empty matrix is returned with the necessary dimensions to store the simulaiton information
#' @export
#'
#' @examples
create_matrix <- function(N_pathogens, N, times, fill=NA){
  m <- lapply(1:N_pathogens, function(x) matrix(fill, nrow=N, ncol=length(times)))
  return(m)
}
