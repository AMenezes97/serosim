#' Calculates An Individual's Current Titre
#'
#' @param boosts_long_infections A matrix containing all long-term boost parameters received from an infection at all possible infection times for each individual and pathogen
#' @param boosts_short_infections A matrix containing all short-term boost parameters received from an infection at all possible infection times for each individual and pathogen
#' @param wanes_long_infections A matrix containing all long-term boost waning parameters from an infection at all possible infection times for each individual and pathogen
#' @param wanes_short_infections A matrix containing all short-term boost waning parameters from an infection at all possible infection times for each individual and pathogen
#' @param boosts_long_vacc A matrix containing all long-term boost parameters received from vaccination at all possible vaccination times for each individual and pathogen
#' @param boosts_short_vacc A matrix containing all short-term boost parameters received from vaccination at all possible vaccination times for each individual and pathogen
#' @param wanes_long_vacc A matrix containing all long-term boost waning parameters from vaccination at all possible vaccination times for each individual and pathogen
#' @param wanes_short_vacc A matrix containing all short-term boost waning parameters from vaccination at all possible vaccination times for each individual and pathogen
#' @param times The total number of time steps for which you want to calculate an individual's titre
#' @param pathogen The pathogen(s) whose titre you want to measure; numeric entries only
#' @param i The individual
#' @param cur_t The current time step
#'
#' @return A titre level gets printed
#' @export
#'
#' @examples
calculate_titre <- function(
  boosts_long_infections, boosts_short_infections,
  wanes_long_infections, wanes_short_infections,
  boosts_long_vacc, boosts_short_vacc,
  wanes_long_vacc, wanes_short_vacc,
  times, pathogen, i, cur_t){

  ## Find gaps of time between now (cur_t) and infection
  t_since_infections <- cur_t - times[which(!is.na(boosts_long_infections[[pathogen]][i,]))]

  ## Extract only finite INFECTION kinetics parameters
  tmp_boosts_long_inf <- boosts_long_infections[[pathogen]][i,]
  tmp_boosts_long_inf <- tmp_boosts_long_inf[!is.na(tmp_boosts_long_inf)]

  tmp_boosts_short_inf <- boosts_short_infections[[pathogen]][i,]
  tmp_boosts_short_inf <- tmp_boosts_short_inf[!is.na(tmp_boosts_short_inf)]

  tmp_wanes_long_inf<- wanes_long_infections[[pathogen]][i,]
  tmp_wanes_long_inf <- tmp_wanes_long_inf[!is.na(tmp_wanes_long_inf)]

  tmp_wanes_short_inf <- wanes_short_infections[[pathogen]][i,]
  tmp_wanes_short_inf <- tmp_wanes_short_inf[!is.na(tmp_wanes_short_inf)]

  ## Add to current titre contribution from infections
  titre <- 0
  for(j in seq_along(tmp_boosts_long_inf)){
    titre <- titre + tmp_boosts_long_inf[j]*max(0, 1-tmp_wanes_long_inf[j]*(t_since_infections[j])) +
      tmp_boosts_short_inf[j]*max(0, 1-tmp_wanes_short_inf[j]*(t_since_infections[j]))
  }

  ## Add to current titre contribution from vaccinations
  t_since_vacc <- cur_t - times[which(!is.na(boosts_long_vacc[[pathogen]][i,]))]

  ## Extract only finite VACCINATION kinetics parameters
  tmp_boosts_long_vacc <- boosts_long_vacc[[pathogen]][i,]
  tmp_boosts_long_vacc <- tmp_boosts_long_vacc[!is.na(tmp_boosts_long_vacc)]

  tmp_boosts_short_vacc <- boosts_short_vacc[[pathogen]][i,]
  tmp_boosts_short_vacc <- tmp_boosts_short_vacc[!is.na(tmp_boosts_short_vacc)]

  tmp_wanes_long_vacc<- wanes_long_vacc[[pathogen]][i,]
  tmp_wanes_long_vacc <- tmp_wanes_long_vacc[!is.na(tmp_wanes_long_vacc)]

  tmp_wanes_short_vacc <- wanes_short_vacc[[pathogen]][i,]
  tmp_wanes_short_vacc <- tmp_wanes_short_vacc[!is.na(tmp_wanes_short_vacc)]

  for(j in seq_along(tmp_boosts_long_vacc)){
    titre <- titre + tmp_boosts_long_vacc[j]*max(0, 1-tmp_wanes_long_vacc[j]*(t_since_vacc[j])) +
      tmp_boosts_short_vacc[j]*max(0, 1-tmp_wanes_short_vacc[j]*(t_since_vacc[j]))
  }
  return(titre)
}
