N <- 10
times <- seq(1,100,by=1)
demography <- tibble(i=1:N,birth=rep(1,N), death=rep(NA,N),location=rep(1,N))
simulation_settings <- list("t_start"=1,"t_end"=max(times))
observation_times <- NULL
lambdas <- array(rep(0.01,length(times)), dim=c(length(times),1,1))
antigen_map <- tibble(exposure_id=1,antigen_id=1)
theta <- list("boost_mean"=2,"boost_sd"=1)

exposure_model <- function(i, t, e, l, lambdas, demography){
    p <- lambdas[t, e, l]
    p
}
immunity_model <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map,...){
    return(1)
}
observation_model <- NULL
draw_parameters <- function(i, t, e, ag, demography, theta, antibody_state, ...){
    boost <- rnorm(1, theta[["boost_mean"]],theta[["boost_sd"]])
    tibble(i=i, t=t, e=e, ag=ag, name="boost",value=boost)    
}
antibody_model <- function(i,t1,ag, exposure_histories,kinetics_parameters,antigen_map){
    exp_history <- exposure_histories[i,1:(t1-1),]
    y <- 0
    if(t1 > 1 & sum(exp_history > 0)){
        tmp_boosts <- kinetics_parameters[[i]] %>% filter(t < t1) %>%
            filter(name == "boost") %>% pull(value)
        for(ts in seq_along(tmp_boosts)){
            y <- y + tmp_boosts[ts]
        }
    }
    y
}

res <- serosim(simulation_settings, demography, observation_times,
        lambdas, antigen_map, theta,
        exposure_model, immunity_model, antibody_model, observation_model, draw_parameters)
image(t(res$antibody_states[,,1]))

