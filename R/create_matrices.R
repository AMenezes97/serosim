#' Title
#'
#' @param N_pathogens
#' @param N
#' @param times
#'
#' @return
#' @export
#'
#' @examples
create_matrices <-function(N_pathogens, N, times){
  infection_histories <- lapply(1:N_pathogens, function(x) matrix(0, nrow=N, ncol=length(times)))
  titres <- lapply(1:N_pathogens, function(x) matrix(0, nrow=N, ncol=length(times)))
  boosts_long_infections <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  boosts_short_infections <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  wanes_long_infections <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  wanes_short_infections <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  boosts_long_vacc <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  boosts_short_vacc <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  wanes_long_vacc <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  wanes_short_vacc <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  all_infection_probs <- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  age<- lapply(1:N_pathogens, function(x) matrix(NA, nrow=N, ncol=length(times)))
  list(infection_histories, titres, boosts_long_infections, boosts_short_infections, wanes_long_infections,
       wanes_short_infections, boosts_long_vacc, boosts_short_vacc, wanes_long_vacc, wanes_short_vacc, all_infection_probs, age)
}

