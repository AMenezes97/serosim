#' Main simulation function for serosim
#' 
#' Simulates a serological survey using custom inputs. The user can specify multiple inputs controlling population demography, simulation timeframe, observation times for each individual, force of infection, and various model functions describing the link between infections and observed antibody titers.
#' 
#' @param simulation_settings
#' @param demography
#' @param observation_times
#' @param lambdas
#' @param antigen_map
#' @param theta
#' @param exposure_model
#' @param immunity_model
#' @param antibody_model
#' @param observation_model
#' @param draw_parameters
#' @param exposure_histories_fixed
#' 
#' @return a list containing the following elements.
#' 
#' @export
#' @examples
serosim <- function(
    ## SIMULATION SETTINGS
    simulation_settings, ## List of parameters governing the simulation settings
    demography=NULL, ## tibble of demographic information for each individual
    observation_times=NULL, ## tibble of observation times and antigen for each individual
    lambdas, ## 3D array giving force of infection for each exposure ID, groups and time
    antigen_map, ## Object determining relationship between exposure IDs and antigens
    theta,
    
    ## FUNCTIONS
    exposure_model, ## calculates the probability of infection given the FOI array, lambdas
    immunity_model, ## function determining probability of infection conditional on lambdas and individuals immune state
    antibody_model, ## function determining antibody state as a function of exposure history and kinetics parameters (theta)
    observation_model, ## function generating observed titers as a function of latent titers and theta
    draw_parameters, ## function to simulate antibody kinetics parameters
    
    ## Pre-specified parameters/events
    exposure_histories_fixed=NULL,
    ...
                    ){
    ## Simulation settings
    t_start <- simulation_settings[["t_start"]]
    t_end <- simulation_settings[["t_end"]]
    times <- seq(simulation_settings[["t_start"]],simulation_settings[["t_end"]],1)
    
    ## Extract key demographic information
    indivs <- unique(demography$i)
    N <- length(indivs)
    
    ## Note "birth" refers to first time point in the population and "removal" refers to time point of removal from population
    if(!("birth" %in% colnames(demography))) {
        demography$birth <- t_start
    } 
    birth_times <- demography %>% select(i, birth) %>% distinct()
    if(!("removal" %in% colnames(demography))) {
        demography$removal <- t_end
    } 
    removal_times <- demography %>% select(i, removal) %>% distinct() 
    
    ## If no groups information provided, assume 1 group
    ## ...
    if(!("group" %in% colnames(demography))) {
        demography$group <- 1
    }
    groups <- demography %>% select(i, group) %>% distinct()
    
    ## Extract information on number of exposure types
    exposure_ids <- unique(antigen_map$exposure_id)
    antigen_ids <- unique(antigen_map$antigen_id)
    N_exposure_ids <- length(exposure_ids)
    N_antigen_ids <- length(antigen_ids)
    
    
    ## Create empty arrays to store exposure histories
    exposure_histories <- array(NA, dim=c(N, length(times), N_exposure_ids))
    exposure_probabilities <- array(NA, dim=c(N, length(times), N_exposure_ids))
    antibody_states <- array(0, dim=c(N, length(times), N_antigen_ids))
    kinetics_parameters <- vector(mode="list",length=N)

    ## Merge in any pre-specified exposure history information
    ## ...
    if(!is.null(exposure_histories_fixed)){
        exposure_histories <- ifelse(!is.na(exposure_histories_fixed), exposure_histories_fixed, exposure_histories) 
    }
    
    message(cat("Beginning simulation\n"))
    ## For each individual
    for(i in indivs){
        message(cat("Individual: ", i, "\n"))
        ## Pull birth time for this individual
        birth_time <- birth_times$birth[i]
        removal_time <- ifelse(is.na(removal_times$removal[i]), simulation_settings[["t_end"]], removal_times$removal[i])
        
        ## Only consider times that the individual was alive for
        simulation_times_tmp <- times[times >= birth_time & times <= removal_time]
        
        ## Go through all times relevant to this individual
        for(t in simulation_times_tmp){
            ## Pull group for this individual at this time 
            g <- as.numeric(demography$group[demography$i==i & demography$times==t]) 
           
             ## Work out antibody state for each antigen
            ## The reason we nest this at the same level as the exposure history generation is
            ## that exposure histories may be conditional on antibody state
            for(ag in antigen_ids){
                antibody_states[i,t,ag] <- antibody_model(i, t, ag, exposure_histories, 
                                                          antibody_states, kinetics_parameters, antigen_map, ...)
            }
            
            ## Work out exposure result for each exposure ID
            for(e in exposure_ids){
                ## Only update if exposure history entry is NA here. If not NA, then pre-specified
                if(is.na(exposure_histories[i,t,e])){
                    ## What is the probability that exposure occurred?
                    prob_exposed <- exposure_model(i, t, e, g, lambdas, demography, ...)
                    
                    ## If an exposure event occurred, what's the probability 
                    ## of successful infection/vaccination?
                    prob_success <- immunity_model(i, t, e, exposure_histories, 
                                                   antibody_states, demography, 
                                                   antigen_map, theta, ...)
                    
                    ## Randomly assign success of exposure event based on immune state
                    successful_exposure <- as.integer(runif(1) < prob_success*prob_exposed)
                    
                    ## Create kinetics parameters for this exposure event
                    ## Each successful exposure event will create a tibble with parameters
                    ## for this event, drawn from information given in theta
                    ## We also pass the demographic information in case we want demography-specific parameters
                    if(successful_exposure == 1){
                        kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                                          draw_parameters(i, t, e, ag, demography, antibody_states, theta, ...))
                    }
                    exposure_histories[i,t,e] <- successful_exposure
                    exposure_probabilities[i,t,e] <- prob_success*prob_exposed
                    if(successful_exposure == 1){
                        for(ag in antigen_ids){ ## Is this correct? Shouldn't it only step through the antigens in the successful exposure? Or is this handled within the function where Ags not in the exposure aren't updated? 
                            antibody_states[i,t,ag] <- antibody_model(i, t, ag, exposure_histories, 
                                                                      antibody_states, kinetics_parameters, antigen_map, ...)
                        }
                    }
                }
            }
        }
    }
    all_kinetics_parameters <- do.call("bind_rows", kinetics_parameters)
    
    ## Reshape antibody states
    antibody_states <- reshape2::melt(antibody_states)
    colnames(antibody_states) <- c("i","t","ag","value")
    antibody_states <- antibody_states %>% arrange(i, t, ag)
    
    ## Reshape exposure histories
    exposure_histories_long <- NULL
    if(sum(exposure_histories, na.rm = TRUE) > 0){
        exposure_histories_long <- reshape2::melt(exposure_histories)
        colnames(exposure_histories_long) <- c("i","t","e","value")
        exposure_histories_long <- exposure_histories_long %>% filter(value != 0) %>% select(-value)
        exposure_histories_long <- exposure_histories_long %>% arrange(i, t, e)
    }
    
    ## Reshape exposure probabilities
    exposure_probabilities_long <- reshape2::melt(exposure_probabilities)
    colnames(exposure_probabilities_long) <- c("i","t","e","value")
    exposure_probabilities_long <- exposure_probabilities_long %>% arrange(i, t, e)
    ## Observation process
    if(!is.null(observation_times)){
        observed_antibody_states <- observation_model(left_join(observation_times,antibody_states), theta, demography, ...)
    } else {
        observed_antibody_states <- observation_model(antibody_states, theta, demography, ...)
    }
    
    return(list("exposure_histories"=exposure_histories,
                "exposure_histories_long"=exposure_histories_long,
                "exposure_probabilities"=exposure_probabilities,
                "exposure_probabilities_long"=exposure_probabilities_long,
                "antibody_states"=antibody_states,
                "observed_antibody_states"=observed_antibody_states,
                "kinetics_parameters"=all_kinetics_parameters))
}
