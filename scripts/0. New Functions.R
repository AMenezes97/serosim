
# generate covariance matrix for points in `x` using given kernel function
cov_matrix <- function(x, kernel_fn, ...) {
  outer(x, x, function(a, b) kernel_fn(a, b, ...))
}

# given x coordinates, take N draws from kernel function at those points
draw_samples <- function(x, N, seed = 1, kernel_fn, ...) {
  Y <- matrix(NA, nrow = length(x), ncol = N)
  #set.seed(seed)
  for (n in 1:N) {
    K <- cov_matrix(x, kernel_fn, ...)
    Y[, n] <- mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
  }
  Y
}

se_kernel <- function(x, y, sigma = 1, length = 1) {
  sigma^2 * exp(- (x - y)^2 / (2 * length^2))
}


titre_protection <- function(titre, alpha1, beta1){
  risk <- 1 - 1/(1 + exp(beta1*(titre - alpha1)))
  risk
}

p_infection <- function(phi, titre, alpha1, beta1){
  p <- phi*(1-titre_protection(titre, alpha1 , beta1))
  p
}


normal_to_lognormal_mean <- function(normmean, normsd) {
  phi <- sqrt(normsd ^ 2 + normmean ^ 2)
  meanlog <- log(normmean ^ 2 / phi)
  return(meanlog)
}

normal_to_lognormal_sd <- function(normmean, normsd) {
  phi <- sqrt(normsd ^ 2 + normmean ^ 2)
  sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
  return(sdlog)
}


draw_parameters <- function(par1, par2, distribution="log-normal"){
  if(par1==0 & par2==0){
    return(0)
  }
  if(distribution=="log-normal"){

    return(rlnorm(1, normal_to_lognormal_mean(par1,par2), normal_to_lognormal_sd(par1,par2)))
  }
  # if(par1==0 & par2==0){
  #   return(0)
  # }
}


simulate_boost <- function(current_titre, boost_mean, boost_sd, titre_ceiling_gradient, titre_ceiling_threshold, distribution){
  boost <- draw_parameters(boost_mean, boost_sd, distribution)
  titre_threshold <- min(current_titre, titre_ceiling_threshold)
  boost <- boost*(1-titre_ceiling_gradient*titre_threshold)
  boost
}


create_infection_exposure_parameters <- function(boosts_long_infections, boosts_short_infections,
                                       wanes_long_infections, wanes_short_infections,
                                       pathogen, i, t, kinetics_pars){
  boosts_long_infections[[pathogen]][i,t] <- simulate_boost(
    current_titre,
    kinetics_pars[kinetics_pars$name=="boost_infection_long_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="boost_infection_long_var" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_gradient" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_threshold" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  boosts_short_infections[[pathogen]][i,t] <- simulate_boost(
    current_titre,
    kinetics_pars[kinetics_pars$name=="boost_infection_short_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="boost_infection_short_var" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_gradient" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="titre_ceiling_threshold" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  wanes_long_infections[[pathogen]][i,t] <- draw_parameters(
    kinetics_pars[kinetics_pars$name=="wane_infection_long_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="wane_infection_long_var" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  wanes_short_infections[[pathogen]][i,t] <- draw_parameters(
    kinetics_pars[kinetics_pars$name=="wane_infection_short_mean" & kinetics_pars$pathogen == pathogen,"value"],
    kinetics_pars[kinetics_pars$name=="wane_infection_short_var" & kinetics_pars$pathogen == pathogen,"value"],
    distribution="log-normal"
  )
  list(boosts_long_infections, boosts_short_infections,
       wanes_long_infections, wanes_short_infections)
}


#this function generates all exposure parameters for short and long term vaccination boosts by calling the simulate_function which samples from log-normal distributions with the specified mean and var
#this is the random effects layer(?)
create_vaccination_exposure_parameters <- function(boosts_long_vacc, boosts_short_vacc,
                                                   wanes_long_vacc, wanes_short_vacc,
                                                   pathogen, i, t, kinetics_pars){
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

# calculate_antibody_trajectory <- function(times, boosts_long, boosts_short, wanes_long, wanes_short){
#   ## Add to current titre contribution from infections
#
#   titres <- numeric(length(times))
#
#   for(i in seq_along(times)){
#     t <- times[i]
#     for(j in seq_along(boosts_long)){
#
#       t_since_infection <- t - times[j]
#         if(t_since_infection >= 0){
#           titres[t] <- titres[t] + boosts_long[j]*max(0, 1-wanes_long[j]*(t_since_infection)) +
#           boosts_short[j]*max(0, 1-wanes_short[j]*(t_since_infection))
#       }
#     }
#   }
#   titres
# }

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
