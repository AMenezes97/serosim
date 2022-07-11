#' Simulate Random Birth Times
#'
#' @param N The number of individuals in the simulation
#' @param times The total number of time steps in the simulation
#' @param limit This number limits the last month an individual is born; if you want to ensure that all individuals simulated are above age of vaccination then enter the time step at which individuals are eligible for vaccination; defaults to 0 which allows individuals to be born up until the second to last step which ensures that each individual is alive for sampling at the last time step
#' @return Random birth times for each individual are returned
#' @export
#'
#' @examples simulate_birth_times(500, 1:100, limit=9) #Simulates random birth times for 500 individuals over 100 time  steps and ensures that all indiivduals are above 9 time steps old by the last time step
simulate_birth_times <- function(N, times, limit=0){
    birth_times <- sample(times[1:(length(times)-(limit+1))], N, replace =TRUE)
    return(birth_times)
}
