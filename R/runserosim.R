#' Main simulation function for serosim
#' 
#' @description Simulates a serological survey using user-specified inputs. The user can specify multiple inputs controlling population demography, simulation time period, observation times for each individual, force of exposure, and various model functions describing the link between infections and observed antibody titers.
#' 
#' @param simulation_settings A list of parameters governing the simulation time step settings
#' @param demography A tibble of relevant demographic information for each individual in the simulation. This tibble only requires 1 column (i) where all individuals in the simulation are listed by row. This is where the sample size for the simulation will be extracted from. If no information is included for birth and removal time, the model will assume that birth time is the initial time point and removal time is the final time point across all individuals. 
#' @param observation_times A tibble of observation times and biomarkers measured for each individual
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time
#' @param biomarker_map A table specifying the relationship between exposure IDs and biomarker IDs
#' @param model_pars A tibble of parameters needed for the antibody kinetics model, immunity model, observation model and the draw_parameters function 
#' @param exposure_model A function which calculates the probability of exposure given the foe_pars array
#' @param immunity_model A function determining the probability of an exposure leading to successful infection or vaccination for a given individual
#' @param antibody_model A function determining the antibody state as a function of infection and vaccination history and antibody kinetics parameters (model_pars)
#' @param observation_model A function generating observed titers as a function of latent titers and kinetics parameters (model_pars)
#' @param draw_parameters A function to simulate parameters antibody kinetics model, immunity model, observation model from model_pars
#' @param exposure_histories_fixed (optional) A 3D array indicating the exposure history (1 = exposed) for each individual (dimension 1) at each time (dimension 2) for each exposure ID (dimension 3).  Here, users can input pre-specified information if exposure histories are known for any individuals.
#' @param VERBOSE (optional) If an integer is specified; an update message will be printed once the simulation reaches that individual and every multiple thereafter; defaults to NULL 
#' 
#' @return a list containing the following elements: force of exposure, exposure probabilities, exposure histories, antibody states, observed antibody states, and kinetics parameters 
#' 
#' @export
#' @examples
runserosim <- function(
    ## SIMULATION SETTINGS
    simulation_settings, ## List of parameters governing the simulation settings
    demography=NULL, ## tibble of demographic information for each individual
    observation_times=NULL, ## tibble of observation times and biomarkers measured for each individual
    foe_pars, ## 3D array giving force of infection for each exposure ID, groups and time
    biomarker_map, ## Object determining relationship between exposure IDs and biomarkers
    model_pars,
    
    ## FUNCTIONS
    exposure_model, ## calculates the probability of infection given the FOI array, foe_pars
    immunity_model, ## function determining probability of infection conditional on foe_pars and individuals immune state
    antibody_model, ## function determining antibody state as a function of exposure history and kinetics parameters (model_pars)
    observation_model, ## function generating observed titers as a function of latent titers and model_pars
    draw_parameters, ## function to simulate antibody kinetics parameters
    
    ## Pre-specified parameters/events (optional)
    exposure_histories_fixed=NULL,
    
    ## UPDATE MESSAGE
    VERBOSE=NULL,
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
    exposure_ids <- unique(biomarker_map$exposure_id)
    biomarker_ids <- unique(biomarker_map$biomarker_id)
    N_exposure_ids <- length(exposure_ids)
    N_biomarker_ids <- length(biomarker_ids)
    
    
    ## Create empty arrays to store exposure histories
    exposure_histories <- array(NA, dim=c(N, length(times), N_exposure_ids))
    exposure_probabilities <- array(NA, dim=c(N, length(times), N_exposure_ids))
    exposure_force <-array(NA, dim=c(N, length(times), N_exposure_ids))
    antibody_states <- array(0, dim=c(N, length(times), N_biomarker_ids))
    kinetics_parameters <- vector(mode="list",length=N)

    ## Merge in any pre-specified exposure history information
    ## ...
    if(!is.null(exposure_histories_fixed)){
        exposure_histories <- ifelse(!is.na(exposure_histories_fixed), exposure_histories_fixed, exposure_histories) 
    }
    
    # message(cat("Beginning simulation\n"))
    ## For each individual
    for(i in indivs){
        ## Print update message
        update(VERBOSE,i)
        ## Pull birth time for this individual
        birth_time <- birth_times$birth[i]
        removal_time <- ifelse(is.na(removal_times$removal[i]), simulation_settings[["t_end"]], removal_times$removal[i])
        
        ## Only consider times that the individual was alive for
        simulation_times_tmp <- times[times >= birth_time & times <= removal_time]
        
        ## Go through all times relevant to this individual
        for(t in simulation_times_tmp){
            ## Pull group for this individual at this time 
            g <- as.numeric(demography$group[demography$i==i & demography$times==t]) 
           
            ## Work out antibody state for each biomarker
            ## The reason we nest this at the same level as the exposure history generation is
            ## that exposure histories may be conditional on antibody state
            for(b in biomarker_ids){
                antibody_states[i,t,b] <- antibody_model(i, t, b, exposure_histories, 
                                                          antibody_states, kinetics_parameters, biomarker_map, ...)
            }
            
            ## Work out exposure result for each exposure ID
            for(x in exposure_ids){
                ## Only update if exposure history entry is NA here. If not NA, then pre-specified
                if(is.na(exposure_histories[i,t,x])){
                    ## What is the probability that exposure occurred?
                    prob_exposed <- exposure_model(i, t, x, g, foe_pars, demography, ...)
                    
                    ## If an exposure event occurred, what's the probability 
                    ## of successful exposure?
                    prob_success <- immunity_model(i, t, x, exposure_histories, 
                                                   antibody_states, demography, 
                                                   biomarker_map, model_pars, ...)

                                        ## Randomly assign success of exposure event based on immune state
                    successful_exposure <- as.integer(runif(1) < prob_success*prob_exposed)
                    
                    ## Simulate kinetics parameters for this exposure event
                    ## Each successful exposure event will create a tibble with parameters
                    ## for this event, drawn from information given in model_pars
                    if(successful_exposure == 1){
                        kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                                          draw_parameters(i, t, x, b, demography, antibody_states, model_pars, ...))
                    }
                    exposure_histories[i,t,x] <- successful_exposure
                    exposure_probabilities[i,t,x] <- prob_success*prob_exposed
                    exposure_force[i,t,x] <- prob_exposed
                    if(successful_exposure == 1){
                        for(b in biomarker_ids){
                            antibody_states[i,t,b] <- antibody_model(i, t, b, exposure_histories, 
                                                                      antibody_states, kinetics_parameters, biomarker_map, ...)
                        }
                    }
                }
            }
        }
    }
    all_kinetics_parameters <- do.call("bind_rows", kinetics_parameters)
    
    ## Reshape antibody states
    antibody_states <- reshape2::melt(antibody_states)
    colnames(antibody_states) <- c("i","t","b","value")
    antibody_states <- antibody_states %>% arrange(i, t, b)
    
    ## Reshape exposure histories
    exposure_histories_long <- NULL
    if(sum(exposure_histories, na.rm = TRUE) > 0){
        exposure_histories_long <- reshape2::melt(exposure_histories)
        colnames(exposure_histories_long) <- c("i","t","x","value")
        # exposure_histories_long <- exposure_histories_long %>% filter(value != 0) %>% select(-value)
        exposure_histories_long <- exposure_histories_long %>% arrange(i, t, x)
    }
    
    ## Reshape probability of a successful exposure event
    exposure_probabilities_long <- reshape2::melt(exposure_probabilities)
    colnames(exposure_probabilities_long) <- c("i","t","x","value")
    exposure_probabilities_long <- exposure_probabilities_long %>% arrange(i, t, x)
    
    ## Reshape exposure probabilities
    exposure_force_long <- reshape2::melt(exposure_force)
    colnames(exposure_force_long) <- c("i","t","x","value")
    exposure_force_long <- exposure_force_long %>% arrange(i, t, x)
    
    ## Observation process
    if(!is.null(observation_times)){
        observed_antibody_states <- observation_model(left_join(observation_times,antibody_states), model_pars, ...)
    } else {
        observed_antibody_states <- observation_model(antibody_states, model_pars, ...)
    }
    
    return(list("exposure_histories"=exposure_histories,
                "exposure_histories_long"=exposure_histories_long,
                "exposure_probabilities"=exposure_probabilities,
                "exposure_probabilities_long"=exposure_probabilities_long,
                "exposure_force"=exposure_force,
                "exposure_force_long"=exposure_force_long,
                "antibody_states"=antibody_states,
                "observed_antibody_states"=observed_antibody_states,
                "kinetics_parameters"=all_kinetics_parameters))
}
