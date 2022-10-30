#' Main simulation function for serosim
#' 
#' Simulates a serological survey using custom inputs. The user can specify multiple inputs controlling population demography, simulation timeframe, observation times for each individual, force of infection, and various model functions describing the link between infections and observed antibody titers.
#' 
#' @param simulation_settings A list of parameters governing the simulation time step settings
#' @param demography A tibble of relevant demographic information for each individual in the simulation. This tibble only requires 1 column (i) where all individuals in the simulation are listed by row. This is where the sample size for the simulation will be extracted from. If no information is included for birth and removal time, the model will assume that birth time is the initial time point and removal time is the final time point across all individuals. 
#' @param observation_times A tibble of observation times and antigen for each individual
#' @param foe_pars A 3D array providing the rate of infection or vaccination for each exposure ID, group and time
#' @param antigen_map An object specifying the relationship between exposure IDs and antigen IDs
#' @param theta A tibble of parameters needed for the antibody kinetics model, immunity model, observation model and the draw_parameters function 
#' @param exposure_model A function which calculates the probability of exposure given the foe_pars array
#' @param immunity_model A function determining the probability of an exposure leading to successful infection or vaccination for a given individual
#' @param antibody_model A function determining the antibody state as a function of infection and vaccination history and antibody kinetics parameters (model_pars)
#' @param observation_model A function generating observed titers as a function of latent titers and kinetics parameters (model_pars)
#' @param draw_parameters A function to simulate parameters antibody kinetics model, immunity model, observation model from model_pars
#' @param exposure_histories_fixed A 3D array indicating the exposure history (1 = exposed) for each individual (dimension 1) at each time (dimension 2) for each exposure ID (dimension 3).  Here, users can input pre-specified information if exposure histories are known for any individuals.
#' 
#' @return a list containing the following elements: exposure probabilities, exposure histories, antibody states, observed antibody states, and kinetics parameters 
#' 
#' @export
#' @examples
runserosim <- function(
    ## SIMULATION SETTINGS
    simulation_settings, ## List of parameters governing the simulation settings
    demography=NULL, ## tibble of demographic information for each individual
    observation_times=NULL, ## tibble of observation times and antigen for each individual
    foe_pars, ## 3D array giving force of infection for each exposure ID, groups and time
    antigen_map, ## Object determining relationship between exposure IDs and antigens
    theta,
    
    ## FUNCTIONS
    exposure_model, ## calculates the probability of infection given the FOI array, foe_pars
    immunity_model, ## function determining probability of infection conditional on foe_pars and individuals immune state
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
    
    # message(cat("Beginning simulation\n"))
    ## For each individual
    for(i in indivs){
        # message(cat("Individual: ", i, "\n"))
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
                    prob_exposed <- exposure_model(i, t, e, g, foe_pars, demography, ...)
                    
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
                        for(ag in antigen_ids){
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
        # exposure_histories_long <- exposure_histories_long %>% filter(value != 0) %>% select(-value)
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
