#' Simulate Vaccine Histories
#'
#' @param N The number of individuals in the simulation
#' @param times The number of time steps in the simulation
#' @param birth_times The time of birth of each individual in the simulation
#' @param vacc_min The age at which an individual becomes eligible for vaccination
#' @param vacc_max The age at which an individual is no longer eligible for vaccination
#' @param prob_vaccination The overall probability of being vaccinated at some point in an individual's life
#' @param vacc_num The maximum number of vaccines an individual can receive; defaults to 1
#'
#' @return A matrix with each individual's vaccination histories is returned
#' @export
#'
#' @examples
simulate_vaccine_histories <- function(N, times, birth_times, vacc_min, vacc_max, prob_vaccination, vacc_num=1){
  vaccine_histories <- matrix(0, nrow=N, ncol=length(times))
  for(i in 1:N){
    tmp <- vaccine_histories[i,]
    #Find vaccination time(s)
    vacc_time <- sample(times[times >= birth_times[i] + (vacc_min) & times <= birth_times[i] + (vacc_max)], vacc_num)
    tmp[vacc_time] <- ifelse(runif(1) < prob_vaccination, 1, 0)
    #Cannot be vaccinated before birth
    tmp[times < birth_times[i]] <- NA
    vaccine_histories[i,] <- tmp
  }
  return(vaccine_histories)
}
