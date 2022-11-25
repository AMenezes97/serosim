#' Generate Population Demography Data Set 
#' 
#' @description Generates a tibble of random demographic information for each individual in the simulation. The returned `demography` tibble can be modified post-hoc to use user-specified distributions and values.
#'
#' @param N The number of individuals in the simulation.
#' @param times A vector of each time step in the simulation.
#' @param birth_times A vector of birth times for each individual; defaults to NULL; if birth_times is not specified then the function will simulate uniformly distributed birth times for each individual from the times vector.
#' @param age_min A number matching the time resolution of `times` giving the youngest age possible by the end of the simulation; defaults to 0 which means individuals can be born up until the penultimate time step.
#' @param removal_min The minimum age at which an individual can be removed from the population. Defaults to 0. 
#' @param removal_max The maximum age at which an individual can be removed from the population. Defaults to max(times).
#' @param prob_removal The probability that an individual will be removed from the population during the simulation, representing e.g., death or study attrition.
#' @param aux An optional list of lists describing additional demographic factors. Each list must contain the variable name, a vector of possible factor levels, and their proportions to simulate from. Note that this is intended for categorical, uncorrelated variables. More complex demographic information should be added post-hoc. Defaults to NULL.
#'
#' @return A tibble of relevant demographic information for each individual in the simulation is returned; this output matches the required `demography` input for the \code{\link{runserosim}} function.
#' @family demography
#' @export
#'
#' @examples 
#' ## Example 1
#' generate_pop_demography(10, times=1:120, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3)
#' 
#' ## Example 2
#' birth_times <- rpois(100, 5)
#' generate_pop_demography(N=100, times=1:120, birth_times=birth_times, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3)
#' 
#' ## Example 3
#' aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "proportion"=c(0.5,0.5)),
#'             "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "proportion"=c(0.25,0.25,0.25,0.25)) )
#' generate_pop_demography(10, 1:120, age_min=0, removal_min=0, removal_max=120, prob_removal=0.3, aux=aux)
generate_pop_demography<-function(N, times, birth_times=NULL, age_min=0, removal_min=min(times), removal_max=max(times), prob_removal, aux=NULL){
    ## If no birth times provided, simulate
    if(is.null(birth_times)){
       birth_times <- simulate_birth_times(N, times, age_min)
    }
    
    ## Simulate removal times
    removal_times <- simulate_removal_times(N, times=times,birth_times=birth_times, removal_min=removal_min, removal_max=removal_max, prob_removal)
    
    demog <- tibble(i=1:N,birth=birth_times,removal=removal_times)
    demog <- demog %>% left_join(expand_grid(i=1:N,time=times))

    ## If there is additional auxiliary demographic information, merge this with default demography table
    if(!is.null(aux)){
        vars <- NULL
        
        for(var in seq_along(aux)){
            ## For each demographic factor, sample a value for each individual
            vars[[var]] <- tibble(sample(aux[[var]][["options"]],N, prob=aux[[var]][["distribution"]],replace=TRUE))
            colnames(vars[[var]])[1] <- aux[[var]]["name"]
        }
        vars <- do.call("bind_cols",vars)
        vars <- vars %>% mutate(i=1:n())
        demog<- demog %>% left_join(vars, by="i")
    }
    return(demog)
}

#' Simulate Random Birth Times
#' 
#' @description This function simulates uniformly distributed birth times for a specified number of individuals from a provided vector
#' 
#' @inheritParams generate_pop_demography
#' @return A vector of simulated birth times for each individual is returned
#' @family demography
#' @export
#'
#' @examples 
#' ## Simulate random birth times for 500 individuals over 100 time steps and ensures that all individuals are above 9 time steps old by the last time step
#' simulate_birth_times(500, 1:100, age_min=9) 
simulate_birth_times <- function(N, times, age_min=0){
    if(age_min >= max(times)){
        message("age_min is greater than or equal to the final time step. Setting to max(times)-1.")
        age_min <- max(times)-1
    }
  birth_times <- sample(times[1:(length(times)-(age_min+1))], N, replace =TRUE)
  return(birth_times)
}


#' Simulate Removal Times for Individuals 
#'
#' @description This function simulates uniformly distributed removal times for a specified number of individuals from the provided times vector. This might represent e.g., death or study attrition.
#' 
#' @inheritParams generate_pop_demography
#' @return A vector of all individual's removal times is returned. NA represents no removal. 
#' @family demography
#' @export
#'
#' @examples
#' ## Simulate random removal times for all individuals; Individuals have a 0.4 probability of being removed at sometime after they are 10 time steps old and before they are 99 time steps old 
#' simulate_removal_times(500,1:100,removal_min=10,removal_max=99, prob_removal=0.4)
simulate_removal_times <- function(N, times, birth_times, removal_min=0, removal_max=max(times), prob_removal=0){
    if(removal_min < min(times)){
        message("removal_min is less than the first time step. Setting to min(times).")
        removal_min <- min(times)
    }
    if(removal_max > max(times)){
        message("removal_max is greater than the final time step. Setting to max(times).")
        removal_max <- max(times)
    }
    
  removal_histories <- matrix(NA, nrow=N, ncol=length(times))
  ## For each individual, randomly sample a removal time with probability prob_removal
  ## between removal_min and removal_max
  sapply(birth_times, function(x) ifelse(runif(1)<prob_removal, sample(max(x, removal_min):removal_max, 1), NA))
}



#' runserosim function update message
#' 
#' @description This function is used within runserosim to produce messages to inform the user which individual the simulation is currently running
#'
#' @param VERBOSE The interval of individuals in which messages will be printed at
#' @param i The individual that the simulation is currently on 
#'
#' @return A printed statement indicating an individual number is returned 
#'
#' @examples
update <- function(VERBOSE, i){ 
  if(!is.null(VERBOSE)){
   if(i %% VERBOSE == 0){
    message(cat("Individual: ", i,"\n"))
  }
  }
}

#' Reformat `biomarker_map` or `model_pars` exposure and biomarker variables
#' 
#' @description This function will reformat the `biomarker_map` or `model_pars` objects so that exposure_ID and biomarker_ID are either both numeric (if passed as characters) or characters (if passed as numeric).
#' 
#' @param input_map A tibble specifying the relationship between exposure IDs and biomarker IDs
#' @param exposure_key Optional vector giving the character-index relationship for `exposure_id`
#' @param biomarker_key Optional vector giving the character-index relationship for `biomarker_id`
#' @return `input_map` is returned with unique numeric or character inputs for exposure and biomarker IDs
#' @export
#'
#' @examples
#' ## Convert characters to numeric
#' biomarker_map <- tibble(exposure_id=c("infection","vaccination"),biomarker_id=c("IgG","IgG"))
#' reformat_biomarker_map(biomarker_map)
#' 
#' ## Convert numeric to characters
#' biomarker_map <- tibble(exposure_id=c(1,2),biomarker_id=c(1,1))
#' reformat_biomarker_map(biomarker_map, exposure_key=c("infection","vaccination"),biomarker_key=c(1))
reformat_biomarker_map<-function(input_map, exposure_key=NULL, biomarker_key=NULL){
    if(any(apply(input_map, 2, class) == "character")){
        ## Converting character to numeric
        input_map$exposure_id <- as.numeric(as.factor(input_map$exposure_id))
        input_map$biomarker_id <- as.numeric(as.factor(input_map$biomarker_id))
    } else {
        ## Converting numeric to character
        if(!is.null(exposure_key) & !is.null(biomarker_key)){
            input_map$exposure_id <- exposure_key[input_map$exposure_id]
            input_map$biomarker_id <- biomarker_key[input_map$biomarker_id]
        }
    }
    input_map
}