antibody_model <- function(i, t, ag, exposure_histories,  antibody_states, kinetics_parameters, antigen_map, ...){
    ## Use exposures for this individuals
    exp_hist <- exposure_histories[i,,]
    ## Get a vector of possible exposure times and types
    times <- 1:nrow(exp_hist)
    exposure_types <- 1:ncol(exp_hist)
    rep_times <- rep(times, ncol(exp_hist))
    rep_exp_type <- rep(exposure_types, length(times))
    
    ## Only consider exposures relevant to this antigen
    use_exposures <- antigen_map[antigen_map$Ag == ag,"E"]
    ## Convert exposure history to a vector with length matching rep_times and rep_exp_type
    exp_hist <- c(exp_hist[,use_exposures])
    
    ## Convert exposure history to a vector of times and IDs
    exp_times <-rep_times[which(exp_hist == 1)]
    exp_types <- rep_exp_type[which(exp_hist == 1)]
    
    ## Only consider exposures which happened prior to this sample
    exp_times <- exp_times[exp_times <= t]
    exp_types <- exp_types[exp_types <= t]
    
    ## Calculate boosts at each infection
    ## If there are infections
    ## Create vector to store effective boost to this antigen from each exposure
    boosts <- numeric(length(exp_times))
    y <- 0
    if(length(exp_times) > 1){
   
        ## Go through each possible exposure
        for(i in seq_along(exp_times)){
            ## Boost without titer mediated boosting -- just unmitigated 
            ## boost with cross reactivity modifier, and considering initial titer
            boosts[i] <- kinetics_pars["boost"] *
                cross_reactivity(exp_types[i], ag, kinetics_pars, antigen_map)
            
            ## If this is not the first boost, then use the previous boosts to calculate titer at 
            ## time of this exposure
            if(i > 1){
                boosts[i] <- boosts[i] * 
                ## Titer from all previous boosts up to this point
                titer_mediation(sum(boosts[1:(i-1)] * pmax((1 - kinetics_pars["wane"] * (exp_times[i] - exp_times[1:(i-1)])),0)), kinetics_pars)
            } 
        }
        ## Latent titer is sum of all previous realized boosts with waning
        y <- boosts * pmax((1- kinetics_pars["wane"] * (t - exp_times)),0)
y <- sum(y)
    }
    
    y
}



n <- 1000
tmax <- 12
nexp <- 3
nag <- 2
rands <- sample(c(0,1),size=n*tmax*nexp,replace=TRUE,prob=c(0.85,0.15))
test_exp_hist <- array(rands, dim = c(n, tmax, nexp))

kinetics_pars <- c("boost"=2,"wane"=0.1)
antigen_map <- data.frame(E=c(1,1,2,2,3,3),Ag=c(1,2,1,2,1,2),par=c(1,0,0,1,0.9,0.75))

cross_reactivity <- function(e, ag, kinetics_pars, antigen_map){
    return(1)
    antigen_map[antigen_map$E == e & antigen_map$Ag == ag,"par"]
}
titer_mediation <- function(y, kinetics_pars){
    1
}

titer_dat <- data.frame(t=rep(rep(c(8,9,10,11,12),each=2),n), ag=rep(c(1,2),n*2),i=rep(1:n,each=2*2))

titer_dat$y <- sapply(1:nrow(titer_dat), function(i) antibody_model(titer_dat$i[i], titer_dat$t[i], titer_dat$ag[i],
                                           test_exp_hist, antibody_states, kinetics_pars, antigen_map))
titer_dat$y1 <- antibody_model_wrapper(titer_dat$i-1, titer_dat$t-1, titer_dat$ag-1,test_exp_hist, kinetics_pars, antigen_map)

titer_dat$diff <- titer_dat$y - titer_dat$y1

microbenchmark::microbenchmark(sapply(1:nrow(titer_dat), function(i) antibody_model(titer_dat$i[i], titer_dat$t[i], titer_dat$ag[i],test_exp_hist, antibody_states, kinetics_pars, antigen_map)),
                               antibody_model_wrapper(titer_dat$i-1, titer_dat$t-1, titer_dat$ag-1,test_exp_hist, kinetics_pars, antigen_map))

microbenchmark::microbenchmark(antibody_model_wrapper(titer_dat$i-1, titer_dat$t-1, titer_dat$ag-1,test_exp_hist, kinetics_pars, antigen_map))

## 3000 -> 2700 -> 400

microbenchmark::microbenchmark(antibody_model_cpp(0, 11, 1, test_exp_hist, kinetics_pars, antigen_map),
                              
antibody_model(1, 12, 1,test_exp_hist, antibody_states, kinetics_pars, antigen_map))



antibody_model_cpp(45-1, 11-1, 1-1, test_exp_hist, kinetics_pars, antigen_map)
antibody_model(45, 11, 1,
               test_exp_hist, antibody_states, kinetics_pars, antigen_map)

y <- array(0, dim=c(n,tmax,nag))
for(i in 1:n){
    for(t in 1:tmax){
        for(ag in 1:nag){
        y[i,t,ag]  <-   antibody_model(i,t,ag,test_exp_hist, antibody_states, kinetics_pars, antigen_map)
        
    }
    }
}
