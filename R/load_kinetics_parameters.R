#' Import And Define Antibody Kinetics Parameters
#'
#' @param boost_infection_long_mean_1 The mean parameter for infection induced long-term boost for pathogen 1; defaults to 0
#' @param boost_infection_long_var_1 The variance parameter for infection induced long-term boost for pathogen 1; defaults to 0
#' @param wane_infection_long_mean_1 The mean parameter for infection induced long-term boost waning for pathogen 1; defaults to 0
#' @param wane_infection_long_var_1 The variance parameter for infection induced long-term boost waning for pathogen 1; defaults to 0
#' @param boost_infection_short_mean_1 The mean parameter for infection induced short-term boost for pathogen 1; defaults to 0
#' @param boost_infection_short_var_1 The variance parameter for infection induced short-term boost for pathogen 1; defaults to 0
#' @param wane_infection_short_mean_1 The mean parameter for infection induced short-term boost waning for pathogen 1; defaults to 0
#' @param wane_infection_short_var_1 The variance parameter for infection induced short-term boost waning for pathogen 1; defaults to 0
#' @param boost_vacc_long_mean_1 The mean parameter for vaccination induced long-term boost for pathogen 1; defaults to 0
#' @param boost_vacc_long_var_1 The variance parameter for vaccination induced long-term boost for pathogen 1; defaults to 0
#' @param wane_vacc_long_mean_1 The mean parameter for vaccination induced long-term boost waning for pathogen 1; defaults to 0
#' @param wane_vacc_long_var_1 The variance parameter for vaccination induced long-term boost waning for pathogen 1; defaults to 0
#' @param boost_vacc_short_mean_1 The mean parameter for vaccination induced short-term boost for pathogen 1; defaults to 0
#' @param boost_vacc_short_var_1 The variance parameter for vaccination induced short-term boost for pathogen 1; defaults to 0
#' @param wane_vacc_short_mean_1 The mean parameter for vaccination induced short-term boost waning for pathogen 1; defaults to 0
#' @param wane_vacc_short_var_1 The variance parameter for vaccination induced short-term boost waning for pathogen 1; defaults to 0
#' @param titre_ceiling_gradient_1 (1-A)/B; Where A is the proportion of the full boost received at or above the titre_ceiling_threshold (B) for pathogen 1; defaults to 0
#' @param titre_ceiling_threshold_1 The maximum titre level an individual can have before their secondary boost is limited in size for pathogen 1; defaults to 0
#' @param boost_infection_long_mean_2 The mean parameter for infection induced long-term boost for pathogen 2; defaults to 0
#' @param boost_infection_long_var_2 The variance parameter for infection induced long-term boost for pathogen 2; defaults to 0
#' @param wane_infection_long_mean_2 The mean parameter for infection induced long-term boost waning for pathogen 2; defaults to 0
#' @param wane_infection_long_var_2 The variance parameter for infection induced long-term boost waning for pathogen 2; defaults to 0
#' @param boost_infection_short_mean_2 The mean parameter for infection induced short-term boost for pathogen 2; defaults to 0
#' @param boost_infection_short_var_2 The variance parameter for infection induced short-term boost for pathogen 2; defaults to 0
#' @param wane_infection_short_mean_2 The mean parameter for infection induced short-term boost waning for pathogen 2; defaults to 0
#' @param wane_infection_short_var_2 The variance parameter for infection induced short-term boost waning for pathogen 2; defaults to 0
#' @param boost_vacc_long_mean_2 The mean parameter for vaccination induced long-term boost for pathogen 2; defaults to 0
#' @param boost_vacc_long_var_2 The variance parameter for vaccination induced long-term boost for pathogen 2; defaults to 0
#' @param wane_vacc_long_mean_2 The mean parameter for vaccination induced long-term boost waning for pathogen 2; defaults to 0
#' @param wane_vacc_long_var_2 The variance parameter for vaccination induced long-term boost waning for pathogen 2; defaults to 0
#' @param boost_vacc_short_mean_2 The mean parameter for vaccination induced short-term boost for pathogen 2; defaults to 0
#' @param boost_vacc_short_var_2 The variance parameter for vaccination induced short-term boost for pathogen 2; defaults to 0
#' @param wane_vacc_short_mean_2 The mean parameter for vaccination induced short-term boost waning for pathogen 2; defaults to 0
#' @param wane_vacc_short_var_2 The variance parameter for vaccination induced short-term boost waning for pathogen 2; defaults to 0
#' @param titre_ceiling_gradient_2 (1-A)/B; Where A is the proportion of the full boost received at or above the titre_ceiling_threshold (B) for pathogen 2; defaults to 0
#' @param titre_ceiling_threshold_2 The maximum titre level an individual can have before their secondary boost is limited in size for pathogen 2; defaults to 0
#'
#' @return A data frame with kinetics parameters for each pathogen is returned
#' @export
#'
#' @examples
load_kinetics_parameters <-function(boost_infection_long_mean_1=0, boost_infection_long_var_1=0,
                                    wane_infection_long_mean_1=0, wane_infection_long_var_1=0,
                                    boost_infection_short_mean_1=0, boost_infection_short_var_1=0,
                                    wane_infection_short_mean_1=0, wane_infection_short_var_1=0,
                                    boost_vacc_long_mean_1=0, boost_vacc_long_var_1=0,
                                    wane_vacc_long_mean_1=0, wane_vacc_long_var_1=0,
                                    boost_vacc_short_mean_1=0,boost_vacc_short_var_1=0,
                                    wane_vacc_short_mean_1=0, wane_vacc_short_var_1=0,
                                    titre_ceiling_gradient_1=0, titre_ceiling_threshold_1=0,
                                    boost_infection_long_mean_2=0,boost_infection_long_var_2=0,
                                    wane_infection_long_mean_2=0, wane_infection_long_var_2=0,
                                    boost_infection_short_mean_2=0, boost_infection_short_var_2=0,
                                    wane_infection_short_mean_2=0, wane_infection_short_var_2=0,
                                    boost_vacc_long_mean_2=0, boost_vacc_long_var_2=0,
                                    wane_vacc_long_mean_2=0,wane_vacc_long_var_2=0,
                                    boost_vacc_short_mean_2=0, boost_vacc_short_var_2=0,
                                    wane_vacc_short_mean_2=0, wane_vacc_short_var_2=0,
                                    titre_ceiling_gradient_2=0, titre_ceiling_threshold_2=0){
  path<-system.file("extdata", "titre_model_parameters.csv", package="serosim")
  kinetics_pars<-readr::read_csv(path, show_col_types = FALSE)
  kinetics_pars<- as.data.frame(kinetics_pars)
  #PATHOGEN 1
  kinetics_pars$value[1] <- boost_infection_long_mean_1
  kinetics_pars$value[2] <- boost_infection_long_var_1
  kinetics_pars$value[3] <- wane_infection_long_mean_1
  kinetics_pars$value[4] <- wane_infection_long_var_1
  kinetics_pars$value[5] <- boost_infection_short_mean_1
  kinetics_pars$value[6] <- boost_infection_short_var_1
  kinetics_pars$value[7] <- wane_infection_short_mean_1
  kinetics_pars$value[8] <- wane_infection_short_var_1
  kinetics_pars$value[9] <- boost_vacc_long_mean_1
  kinetics_pars$value[10] <- boost_vacc_long_var_1
  kinetics_pars$value[11] <- wane_vacc_long_mean_1
  kinetics_pars$value[12] <- wane_vacc_long_var_1
  kinetics_pars$value[13] <- boost_vacc_short_mean_1
  kinetics_pars$value[14] <- boost_vacc_short_var_1
  kinetics_pars$value[15] <- wane_vacc_short_mean_1
  kinetics_pars$value[16] <- wane_vacc_short_var_1
  kinetics_pars$value[17] <- titre_ceiling_gradient_1
  kinetics_pars$value[18] <- titre_ceiling_threshold_1

  #PATHOGEN 2
  kinetics_pars$value[19] <- boost_infection_long_mean_2
  kinetics_pars$value[20] <- boost_infection_long_var_2
  kinetics_pars$value[21] <- wane_infection_long_mean_2
  kinetics_pars$value[22] <- wane_infection_long_var_2
  kinetics_pars$value[23] <- boost_infection_short_mean_2
  kinetics_pars$value[24] <- boost_infection_short_var_2
  kinetics_pars$value[25] <- wane_infection_short_mean_2
  kinetics_pars$value[26] <- wane_infection_short_var_2
  kinetics_pars$value[27] <- boost_vacc_long_mean_2
  kinetics_pars$value[28] <- boost_vacc_long_var_2
  kinetics_pars$value[29] <- wane_vacc_long_mean_2
  kinetics_pars$value[30] <- wane_vacc_long_var_2
  kinetics_pars$value[31] <- boost_vacc_short_mean_2
  kinetics_pars$value[32] <- boost_vacc_short_var_2
  kinetics_pars$value[33] <- wane_vacc_short_mean_2
  kinetics_pars$value[34] <- wane_vacc_short_var_2
  kinetics_pars$value[35] <- titre_ceiling_gradient_2
  kinetics_pars$value[36] <- titre_ceiling_threshold_2
  return(kinetics_pars)
}
