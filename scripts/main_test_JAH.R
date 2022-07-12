N <- 100
N_exposure_ids <- 2
times <- seq(1,12,by=1)
demography <- tibble(i=1:N,birth=rep(1,N), death=rep(NA,N),location=rep(1,N))
simulation_settings <- list("t_start"=1,"t_end"=max(times))
observation_times <- tibble(i=rep(1:N,each=3), t=sample(8:12, size=N*3,replace=TRUE))
lambdas <- array(c(rep(0.1,length(times)),rep(0.05,length(times))), dim=c(length(times),N_exposure_ids,1))

antigen_map <- tibble(exposure_id=1:2,antigen_id=1:2)

theta <- list("boost_mean"=2,"boost_sd"=1,"obs_sd"=0.25)




exposure_model <- function(i, t, e, g, lambdas, demography){
    p <- 1-exp(-lambdas[t, e, g])
    p
}
immunity_model <- function(i, t, e, exposure_histories, 
                           antibody_states, demography, antigen_map,...){
    return(1)
}

draw_parameters <- function(i, t, e, ag, demography, antibody_states, theta, ...){
    #boost <- rnorm(1, theta[["boost_mean"]],theta[["boost_sd"]])
    boost <- theta[["boost_mean"]]
    tibble(i=i, t=t, e=e,ag=ag,name="boost",value=boost)    
}
antibody_model <- function(i,t1,ag, exposure_histories,antibody_states, kinetics_parameters,antigen_map){
    #return(0)
    ## Only look at exposures relevant to this antigen
    use_exposures <- antigen_map %>% filter(antigen_id == ag) %>% pull(exposure_id)
    exp_history <- exposure_histories[i,1:(t1-1),use_exposures]
    y <- 0
    if(t1 > 1 & sum(exp_history > 0)){
        ## Go through each exposure ID and add its contribution to this antigen
        
        tmp_boosts <- kinetics_parameters[[i]] %>% filter(t < t1) %>%
            filter(name == "boost") %>% pull(value)
        for(ts in seq_along(tmp_boosts)){
            y <- y + tmp_boosts[ts]
        }
    }
    y
}

observation_model <- function(antibody_states, theta, demography){
    antibody_states$observed <- rnorm(nrow(antibody_states),antibody_states$value,theta[["obs_sd"]])
    antibody_states
}


exposure_histories_fixed <- array(NA, dim=c(N, length(times), N_exposure_ids))
exposure_histories_fixed[1:3,1:6,] <- 0

Rprof(tmp<-tempfile())
res <- serosim(simulation_settings, demography, observation_times,
        lambdas, antigen_map, theta,
        exposure_model, immunity_model, antibody_model, observation_model, draw_parameters,
        exposure_histories_fixed=exposure_histories_fixed)
Rprof(NULL)
summaryRprof(tmp)
ggplot(res$antibody_states) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~ag)
ggplot(res$exposure_probabilities_long) + geom_tile(aes(x=t,y=i,fill=value)) + facet_wrap(~e)
ggplot(res$observed_antibody_states) + geom_jitter(aes(x=t,y=value),height=0.1,width=0.25) + facet_wrap(~ag) + scale_x_continuous(limits=range(times))



