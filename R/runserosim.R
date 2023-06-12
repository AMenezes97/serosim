#' Main simulation function for _serosim_
#' 
#' @description Simulates a serological survey using user-specified inputs and models.
#' 
#' @details The `runserosim` function is designed as a plug-and-play multi-level model representing the different stages of the data-generating process for serological data. The user can specify inputs controlling population demography, simulation time period, observation times for each individual, force of exposure, and various model functions describing the link between infections and observed biomarker quantities. See README for full details of the workflow.
#' 
#' There are two optional flags to attempt speed ups with serosim: `attempt_precomputation` and `parallel`. `attempt_precomputation` will check how many interchangeable individuals there are and pre-compute exposure probabilities for these individuals if that would be quicker than solving for each individual separately. `parallel` uses the `foreach` package with `dopar` to run _serosim_ for blocks of individuals on a socket cluster. The number of cores is set with the `n_cores` argument. If any additional packages are needed for the model, then the `parallel_packages` argument should be set.
#' 
#' @param simulation_settings A list of parameters governing the simulation time step settings. Should contain integer entries for _t_start_ and _t_end_
#' @param demography A tibble of demographic information for each individual in the simulation. At a minimum this tibble requires 1 column (i) where all individuals in the simulation are listed by row. This is used to calculate the sample population size. Additional variables can be added by the user, e.g., birth and removal times, see \code{\link{generate_pop_demography}} If not specified, the model will assume that birth time is the initial time point and removal time is the final time point across all individuals
#' @param observation_times A tibble for observation times with three variables: 1) i: individual; 2) t: the timepoint of an observation; 3) b: the biomarker being measured. Defaults to NULL, in which case all latent biomarker states are returned
#' @param foe_pars Any object class (usually a 3D array, but may be a list or other object) providing the force of exposure for each exposure ID, group and time, or parameters of the model generating the FOEs
#' @param biomarker_map A tibble specifying the relationship between exposure IDs and biomarker IDs with two variables: 1) exposure_id: the numeric exposure ID; 2) biomarker_id: the numeric biomarker ID
#' @param model_pars A tibble of parameters needed for the antibody kinetics model, immunity model, observation model and the draw_parameters function 
#' @param exposure_model A function calculating the probability of exposure given the foe_pars array
#' @param immunity_model A function determining the probability of an exposure leading to successful infection or vaccination for a given individual
#' @param antibody_model A function determining the latent biomarker kinetics as a function of immune history and within-host kinetics parameters (\code{model_pars})
#' @param observation_model A function generating observed biomarkers as a function of latent biomarkers and observation model parameters (\code{model_pars})
#' @param draw_parameters A function to simulate parameters for the antibody model, immunity model and observation model from \code{model_pars}
#' @param immune_histories_fixed (optional) Defaults to NULL. Otherwise a 3D array indicating the immune history (1 = exposed, 0 = not exposed) for each individual (dimension 1) at each time (dimension 2) for each exposure ID (dimension 3). Here, users can input pre-specified information if immune histories are known for an individual
#' @param VERBOSE (optional) Defaults to NULL. An integer specifying the frequency at which simulation progress updates are printed
#' @param attempt_precomputation If TRUE, attempts to perform as much pre-computation as possible for the exposure model to speed up the main simulation code. If FALSE, skips this step, which is advised when the simulation is expected to be very fast anyway
#' @param parallel If TRUE, attempts to run _serosim_ using parallel processes with the foreach and parallel packages. Defaults to FALSE
#' @param n_cores The number of cores to use if running in parallel. Defaults to 4
#' @param parallel_packages a character vector of R packages required for the models. Only needed if additional non-imported packages are required.
#' @param ... Any additional arguments needed for the models
#' 
#' @return a list containing the following elements: force of exposure, exposure probabilities, immune histories, antibody states, observed antibody states, and kinetics parameters 
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange 
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom abind abind
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach '%dopar%'
#' @importFrom foreach foreach
#' @importFrom parallel stopCluster
#' @export
#' @examples 
#' ## See package README.
runserosim <- function(
    ## SIMULATION SETTINGS
    simulation_settings, ## List of parameters governing the simulation settings
    demography=NULL, ## tibble of demographic information for each individual
    observation_times=NULL, ## tibble of observation times and biomarkers measured for each individual
    foe_pars, ## any object (usually a 3D array) giving force of infection for each exposure ID, groups and time
    biomarker_map, ## Object determining relationship between exposure IDs and biomarkers
    model_pars,
    
    ## FUNCTIONS
    exposure_model, ## calculates the probability of infection given the FOE array, foe_pars
    immunity_model, ## function determining probability of infection conditional on foe_pars and individual's immune state
    antibody_model, ## function determining biomarker quantity, commonly antibody state, as a function of exposure history and kinetics parameters (model_pars)
    observation_model, ## function generating observed biomarker quantity as a function of latent biomarker quantities and model_pars
    draw_parameters, ## function to simulate biomarker (e.g, antibody)  kinetics parameters
    
    ## Pre-specified parameters/events (optional)
    immune_histories_fixed=NULL,
    
    ## UPDATE MESSAGE
    VERBOSE=NULL,
    attempt_precomputation=TRUE,
    parallel=FALSE,
    n_cores=4,
    parallel_packages=NULL,
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
    groups <- demography %>% select(group) %>% distinct() %>% pull(group)
    N_groups <- length(groups)
    
    ## Make sure demography is arranged by i and times (if there is a time element). 
    ## This is important for subsetting below. 
    ## Also pre-compute which group the individual is in for each time point to speed up
    ## indexing below
    if("times" %in% colnames(demography)){
        use_time <- TRUE
        demography <- demography %>% arrange(i, times)
        ## Create groups matrix
        all_groups <- matrix(NA, nrow=N,ncol=length(times))
        indices <- convert_indices_matrix_to_vector(demography$i, demography$times, N)
        all_groups[indices] <- demography$group
        if(any(is.na(all_groups))) message(cat("Warning - could not find group membership for all individuals at all times. This may lead to an error or unexpected output"))
    } else {
        use_time <- FALSE
        demography <- demography %>% arrange(i)
        all_groups <- demography$group
        all_groups <- matrix(rep(all_groups,each=length(times)),nrow=N,ncol=length(times),byrow=TRUE)
    }
    ## Extract information on number of exposure types
    exposure_ids <- unique(biomarker_map$exposure_id)
    biomarker_ids <- unique(biomarker_map$biomarker_id)
    N_exposure_ids <- length(exposure_ids)
    N_biomarker_ids <- length(biomarker_ids)

    ########################################################################
    ## Checking if pre-compuation of exposure probabilities is possible
    if(!is.null(VERBOSE)) message(cat("Checking for possible pre-computation to save time...\n"))
    if(attempt_precomputation){
        precomputations <- precomputation_checks(N, times,exposure_ids, groups,
                                                 exposure_model, foe_pars, demography, 
                                                 VERBOSE, ...)
        if(precomputations$flag == TRUE){
            foe_pars <- precomputations$foe
            exposure_model <- exposure_model_indiv_fixed
        } 
        if(!is.null(VERBOSE)) message(cat("Using pre-computed exposure probabilities\n"))
        
    } 
    ## If successful precomputation, change the exposure model to just read directly from foe_pars
    ########################################################################
    
    ## Parallelize here
    serosim_internal <- function(tmp_indivs,...){
      ## Create empty arrays to store immune histories
      immune_histories <- array(NA, dim=c(N, length(times), N_exposure_ids))
      exposure_probabilities <- array(NA, dim=c(N, length(times), N_exposure_ids))
      exposure_force <-array(NA, dim=c(N, length(times), N_exposure_ids))
      biomarker_states <- array(0, dim=c(N, length(times), N_biomarker_ids))
      kinetics_parameters <- vector(mode="list",length=N)
      
      ## Merge in any pre-specified immune history information
      ## ...
      if(!is.null(immune_histories_fixed)){
          immune_histories <- ifelse(!is.na(immune_histories_fixed), immune_histories_fixed, immune_histories) 
      }
      
      if(!is.null(VERBOSE)) message(cat("Beginning simulation\n"))
  
      ## For each individual
      for(i in tmp_indivs){
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
              g <- all_groups[i, t]
  
              ## Work out antibody state for each biomarker
              ## The reason we nest this at the same level as the immune history generation is
              ## that immune histories may be conditional on antibody state
              for(b in biomarker_ids){
                  biomarker_states[i,t,b] <- antibody_model(i, t, b, immune_histories, 
                                                            biomarker_states, kinetics_parameters, biomarker_map, ...)
              }
              
              ## Work out exposure result for each exposure ID
              for(x in exposure_ids){
                  ## Only update if immune history entry is NA here. If not NA, then pre-specified
                  if(is.na(immune_histories[i,t,x])){
                      ## What is the probability that exposure occurred?
                      prob_exposed <- exposure_model(i, t, x, g, foe_pars, demography, ...)
                      
                      ## If an exposure event occurred, what's the probability 
                      ## of successful exposure?
                      prob_success <- immunity_model(i, t, x, immune_histories, 
                                                     biomarker_states, demography, 
                                                     biomarker_map, model_pars, ...)
  
                      ## Randomly assign success of exposure event based on immune state
                      successful_exposure <- as.integer(runif(1) < prob_success*prob_exposed)
                      
                      ## Simulate kinetics parameters for this exposure event
                      ## Each successful exposure event will create a tibble with parameters
                      ## for this event, drawn from information given in model_pars
                      if(successful_exposure == 1){
                          kinetics_parameters[[i]] <- bind_rows(kinetics_parameters[[i]],
                                                            draw_parameters(i, t, x, demography, biomarker_states, model_pars, ...))
                      }
                      immune_histories[i,t,x] <- successful_exposure
                      exposure_probabilities[i,t,x] <- prob_success*prob_exposed
                      exposure_force[i,t,x] <- prob_exposed
                      if(successful_exposure == 1){
                          for(b in biomarker_ids){
                              biomarker_states[i,t,b] <- antibody_model(i, t, b, immune_histories, 
                                                                        biomarker_states, kinetics_parameters, biomarker_map, ...)
                          }
                      }
                  }
              }
          }
      }
      return(list(array(biomarker_states[tmp_indivs,,],dim=c(length(tmp_indivs), length(times),N_biomarker_ids)), 
                  kinetics_parameters[tmp_indivs], 
                  immune_histories[tmp_indivs,,], 
                  exposure_probabilities[tmp_indivs,,], 
                  exposure_force[tmp_indivs,,]))
    }
    ## Run either the entire simulation in one go, or split into n_cores jobs and run in parallel
    if(!parallel){
      res <- serosim_internal(indivs,...)
      biomarker_states <- res[[1]]
      kinetics_parameters <- res[[2]]
      immune_histories <- res[[3]]
      exposure_probabilities <- res[[4]]
      exposure_force <- res[[5]]
    } else {
      ## Set up socket cluster to split the population across n_cores processes
      cluster <- makeCluster(n_cores) 
      registerDoParallel(cluster)
      
      if(!is.null(VERBOSE)) message(cat("Running jobs in parallel across ", n_cores, " cores\n"))
      
      n_indivs_per_block <- ceiling(length(indivs)/n_cores)
      indiv_blocks <- split(indivs, ceiling(seq_along(indivs)/n_indivs_per_block))
      
      ## Submit job. Note this is likely where any object or package exports may become an issue.
      res <- foreach(block = 1:n_cores, .packages=c("abind",parallel_packages)) %dopar% {
        serosim_internal(indiv_blocks[[block]],...)
      }
      
      ## Close socket cluster
      stopCluster(cluster)
      
      ## Combine outputs from each process. This could probably be improved without requiring abind.
      biomarker_states <- do.call("abind",args=list(lapply(res, function(x) x[[1]]), along=1))
      immune_histories <- do.call("abind",args=list(lapply(res, function(x) x[[3]]), along=1))
      exposure_probabilities <- do.call("abind",args=list(lapply(res, function(x) x[[4]]), along=1))
      exposure_force <- do.call("abind",args=list(lapply(res, function(x) x[[5]]), along=1))
      kinetics_parameters <- do.call("bind_rows",lapply(res, function(x) x[[2]]))
    }
    
    if(!is.null(VERBOSE)) message(cat("Simulation complete! Cleaning up...\n"))
    
    all_kinetics_parameters <- do.call("bind_rows", kinetics_parameters)
    
    ## Reshape antibody states
    biomarker_states <- reshape2::melt(biomarker_states)
    colnames(biomarker_states) <- c("i","t","b","value")
    biomarker_states <- biomarker_states %>% arrange(i, t, b)
    
    ## Reshape immune histories
    immune_histories_long <- NULL
    if(sum(immune_histories, na.rm = TRUE) > 0){
        immune_histories_long <- reshape2::melt(immune_histories)
        colnames(immune_histories_long) <- c("i","t","x","value")
        # immune_histories_long <- immune_histories_long %>% filter(value != 0) %>% select(-value)
        immune_histories_long <- immune_histories_long %>% arrange(i, t, x)
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
        observed_biomarker_states <- observation_model(left_join(observation_times,biomarker_states), model_pars, ...)
    } else {
        observed_biomarker_states <- observation_model(biomarker_states, model_pars, ...)
    }
    
    ## Remove latent states before individual was born and after they left the study
    biomarker_states <- biomarker_states %>% 
      left_join(demography %>% select(c(i, birth,removal)) %>% distinct()) %>%
      mutate(value=ifelse(t >= birth & t <= removal, value, NA))%>% 
      select(-c(birth,removal))

    ## Remove observations before individual was born and after they left the study
    observed_biomarker_states <- observed_biomarker_states %>% 
      left_join(demography %>% select(c(i, birth,removal)) %>% distinct()) %>%
      mutate(value=ifelse(t >= birth & t <= removal, value, NA))%>% 
      mutate(observed=ifelse(t >= birth & t <= removal, observed, NA))%>% 
      select(-c(birth,removal))
    
    return(list("immune_histories"=immune_histories,
                "immune_histories_long"=immune_histories_long,
                "exposure_probabilities"=exposure_probabilities,
                "exposure_probabilities_long"=exposure_probabilities_long,
                "exposure_force"=exposure_force,
                "exposure_force_long"=exposure_force_long,
                "biomarker_states"=biomarker_states,
                "observed_biomarker_states"=observed_biomarker_states,
                "kinetics_parameters"=all_kinetics_parameters,
                "demography"=demography))
}
