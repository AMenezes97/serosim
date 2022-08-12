#' Find Number Of Individuals Alive At Each Time
#'
#' @param N The number of individuals in the simulation
#' @param times The total number of time steps in the simulation
#'
#' @return A vector with the number of individuals alive at each time is returned
#' @export
#'
#' @examples
find_N_alive <- function(N, times){
  birth_times <- sample(times[1:(length(times)-1)], N, replace =TRUE)
  N_alive <- sapply(times, function (x) length(birth_times[birth_times <= x]))
  return(N_alive)
}