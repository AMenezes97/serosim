#' Create Antibody Kinetics Parameters For Vaccination Boosts
#'
#' @param boosts_long_vacc  Empty matrices for all possible vaccination times for each pathogen and individual to store infection long-term boost parameters
#' @param boosts_short_vacc Empty matrices for all possible vaccination times for each pathogen and individual to store infection short-term boost parameters
#' @param wanes_long_vacc  Empty matrices for all possible vaccination times for each pathogen and individual to store infection long-term boost waning parameters
#' @param wanes_short_vacc Empty matrices for all possible vaccination times for each pathogen and individual to store infection short-term boost wanning parameters
#' @param pathogen The pathogen(s) whose infection parameters you want to generate; numeric entries only
#' @param i Individual
#' @param t The current time step
#' @param kinetics_pars A file with the antibody kinetics parameters
#' @param current_titre An individual's current titre
#'
#' @return The 4 matrices are returned with parameters for all possible vaccination times for each individual and pathogen
#' @export
#'
#' @examples
create_vaccination_exposure_parameters <- function(boosts_long_vacc, boosts_short_vacc,
                                                   wanes_long_vacc, wanes_short_vacc,
                                                   pathogen, i, t, kinetics_pars, current_titre){
  boosts_long_vacc[[pathogen]][i,t] <- simulate_boost(
    current_titre,
    kinetics_pars[kinetics_pars$name=="boost_vacc_long_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="boost_vacc_long_var" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_gradient" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_threshold" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  boosts_short_vacc[[pathogen]][i,t] <- simulate_boost(
    current_titre,
    kinetics_pars[kinetics_pars$name=="boost_vacc_short_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="boost_vacc_short_var" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_gradient" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_threshold" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  wanes_long_vacc[[pathogen]][i,t] <- draw_parameters(
    kinetics_pars[kinetics_pars$name=="wane_vacc_long_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="wane_vacc_long_var" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  wanes_short_vacc[[pathogen]][i,t] <- draw_parameters(
    kinetics_pars[kinetics_pars$name=="wane_vacc_short_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="wane_vacc_short_var" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  list(boosts_long_vacc, boosts_short_vacc,
       wanes_long_vacc, wanes_short_vacc)
}
