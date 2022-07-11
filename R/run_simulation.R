#' Run Simulation To Generate Infection And Vaccination Events
#'
#' @param birth_times Each individual's time of birth
#' @param N The number of individuals in the simulation
#' @param times The number of time steps in the simulation
#' @param N_pathogens The number of pathogens in the simulation
#' @param kinetics_pars A file with the antibody kinetics parameters
#' @param vaccine_histories A matrix with individual's vaccine history at each time step
#' @param pathogen Pathogen
#' @param i Individual
#' @param t Time step
#'
#' @return 13 completed matrices are returned
#' @export
#'
#' @examples
run_simulation<-function(birth_times, N, times, N_pathogens, kinetics_pars, vaccine_histories, pathogen, i, t){

  infection_histories <- create_matrix(N_pathogens, N, times, fill=0)
  titres <- create_matrix(N_pathogens, N, times, fill=0)
  boosts_long_infections <- create_matrix(N_pathogens, N, times)
  boosts_short_infections <- create_matrix(N_pathogens, N, times)
  wanes_long_infections <- create_matrix(N_pathogens, N, times)
  wanes_short_infections <- create_matrix(N_pathogens, N, times)
  boosts_long_vacc <- create_matrix(N_pathogens, N, times)
  boosts_short_vacc <- create_matrix(N_pathogens, N, times)
  wanes_long_vacc <- create_matrix(N_pathogens, N, times)
  wanes_short_vacc <- create_matrix(N_pathogens, N, times)
  all_infection_probs <- create_matrix(N_pathogens, N, times)
  age <- create_matrix(N_pathogens, N, times)

for(i in 1:N){ ## Go through individual
  birth_time <- birth_times[i]
  for(t in times){ ## Step through time
    for(pathogen in 1:N_pathogens){
      ## If not born yet, no titres possible and cannot be infected
      if(t < birth_time){
        infection_histories[[pathogen]][i,t] <- NA
        titres[[pathogen]][i,t] <- NA
        age[[pathogen]][i,t] <- NA
      } else {
        vaccine_history_tmp <- vaccine_histories[i,]
        age_tmp<- t-birth_time
        age[[pathogen]][i,t] <- age_tmp
        ## Get current titre, which will dictate infection probability and titre ceiling effects
        ##calculate titre function adds natural infection and vaccination
        current_titre <- calculate_titre(boosts_long_infections, boosts_short_infections,
                                         wanes_long_infections, wanes_short_infections,
                                         boosts_long_vacc, boosts_short_vacc,
                                         wanes_long_vacc, wanes_short_vacc,
                                         times, pathogen, i, t)

        ## Generate vaccination boosting exposure parameters
        if(vaccine_history_tmp[t] == 1){
          tmp <- create_vaccination_exposure_parameters(boosts_long_vacc, boosts_short_vacc,
                                                        wanes_long_vacc, wanes_short_vacc,
                                                        pathogen, i, t, kinetics_pars, current_titre)
          boosts_long_vacc <- tmp[[1]]
          boosts_short_vacc <- tmp[[2]]
          wanes_long_vacc <- tmp[[3]]
          wanes_short_vacc <- tmp[[4]]
        }


        # Find probability of infection given current titre level and protective effect of titres
        prob_infection <- p_infection(prob_infections[t, pathogen], current_titre, titre_prot_midpoint[pathogen], titre_prot_width[pathogen])
        all_infection_probs[[pathogen]][i, t] <- prob_infection
        infected <- as.integer(runif(1) < prob_infection)


        ## Simulate kinetics parameters from this exposure
        if(infected){
          tmp <- create_infection_exposure_parameters(boosts_long_infections, boosts_short_infections,
                                                      wanes_long_infections, wanes_short_infections,
                                                      pathogen, i, t, kinetics_pars, current_titre)
          boosts_long_infections <- tmp[[1]]
          boosts_short_infections <- tmp[[2]]
          wanes_long_infections <- tmp[[3]]
          wanes_short_infections <- tmp[[4]]
        }

        infection_histories[[pathogen]][i,t] <- infected
        titres[[pathogen]][i,t] <- current_titre
      }
    }
  }
}
  return(list("infection_histories" =infection_histories,  "titres"=titres, "age"=age,
              "vaccine_histories"=vaccine_histories, "boosts_long_infections"=boosts_long_infections,
              "boosts_short_infections"= boosts_short_infections,  "wanes_long_infections"= wanes_long_infections,
              "wanes_short_infections"=wanes_short_infections, "boosts_long_vacc"=boosts_long_vacc,
              "boosts_short_vacc"=boosts_short_vacc, "wanes_long_vacc"=wanes_long_vacc,
              "wanes_short_vacc"=wanes_short_vacc, "all_infection_probs"=all_infection_probs))
}
