#' Run Full Simulation To Generate Infection And Vaccination Events
#'
#' @param N The number of individuals in the simulation
#' @param years The number of years to run the simulation for
#' @param t_periods_per_year The number of time periods per year
#' @param N_pathogens The number of pathogens in the simulation
#' @param obs_time The time step when the serosurvey is conducted
#' @param vacc_min The time step when an individual becomes eligible for vaccination
#' @param vacc_max The time step at which an individual becomes too old to get vaccinated
#' @param kinetics_pars A file with the antibody kinetics parameters
#' @param overall_infection_probs The overall infection probability for each pathogen
#' @param titre_prot_midpoint The titre value at which you have 50% protection from infection for each pathogen
#' @param titre_prot_width A parameter forthe titre protection curve
#' @param prob_vaccination The overall probability of vaccination
#' @param pathogen Pathogen
#' @param i Individual
#' @param t Time step
#'
#' @return
#' @export
#'
#' @examples
run_full_simulation<-function(N, years, t_periods_per_year, N_pathogens, obs_time, vacc_min, vacc_max, kinetics_pars,
                              overall_infection_probs, titre_prot_midpoint, titre_prot_width, prob_vaccination,
                              pathogen, i, t){


  times<- seq(1,years*t_periods_per_year, by=1)
  birth_times <-simulate_birth_times(N, times, limit=vacc_min)
  N_alive <- find_N_alive(N, times)

  #Option 1: Random force of infection for all pathogens
  FOIs <- draw_samples(times, N_pathogens, kernel_fn = se_kernel, sigma=2,length = 24)
  # ?generate_prob_infections
  prob_infections<-generate_prob_infections(FOIs, N_pathogens, overall_infection_probs)

  #Simulate vaccine histories
  # ?simulate_vaccine_histories
  vaccine_histories<-simulate_vaccine_histories(N, times, birth_times, vacc_min, vacc_max, prob_vaccination)

  #Create empty vectors to store simulation outputs
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
  return(list("times"=times, "birth_times"=birth_times, "N_alive"=N_alive, "FOIs"=FOIs, "prob_infections"=prob_infections,
              "infection_histories" =infection_histories,  "titres"=titres, "age"=age,
              "vaccine_histories"=vaccine_histories, "boosts_long_infections"=boosts_long_infections,
              "boosts_short_infections"= boosts_short_infections,  "wanes_long_infections"= wanes_long_infections,
              "wanes_short_infections"=wanes_short_infections, "boosts_long_vacc"=boosts_long_vacc,
              "boosts_short_vacc"=boosts_short_vacc, "wanes_long_vacc"=wanes_long_vacc,
              "wanes_short_vacc"=wanes_short_vacc, "all_infection_probs"=all_infection_probs))
}
