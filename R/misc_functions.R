#' Simulate Random Birth Times
#' 
#' @description This function simulates random birth times for a specified number of individuals across a range of times
#' 
#' @param N The number of individuals in the simulation
#' @param times  A vector of each time step in the simulation
#' @param limit A number indicating the youngest age possible by the end of the simulation; defaults to 0 which means individuals can be born up until the second to last time step
#'
#' @return A vector of simulated birth times for each individual is returned
#' @export
#'
#' @examples 
#' ## Simulate random birth times for 500 individuals over 100 time steps and ensures that all individuals are above 9 time steps old by the last time step
#' simulate_birth_times(500, 1:100, limit=9) 
simulate_birth_times <- function(N, times, limit=0){
  birth_times <- sample(times[1:(length(times)-(limit+1))], N, replace =TRUE)
  return(birth_times)
}


#' Simulate  Removal Times for Individuals 
#'
#' @description This function simulates random removal times for a specified number of individuals across a range of times
#' 
#' @param N The number of individuals in the simulation
#' @param times  A vector of each time step in the simulation
#' @param birth_times A vector of all individual's birth times
#' @param removal_min The minimum age at which an individual can be removed from the population 
#' @param removal_max The maximum age at which an individual can be removed from the population 
#' @param prob_removal The probability that an individual will be removed from the population
#'
#' @return A vector of all individual's removal times is returned 
#' @export
#'
#' @examples
#' ## First, simulate random birth times for 500 individuals over 100 time steps and ensures that all individuals are above 9 time steps old by the last time step
#' births<-simulate_birth_times(500, 1:100, limit=9) 
#' ## Simulate random removal times for all individuals; Individuals have a 0.4 probability of being removed at sometime after they are 10 time steps old and before they are 99 time steps old 
#' simulate_removal_times(500,1:100,birth_times=births,removal_min=10,removal_max=99, prob_removal=0.4)
simulate_removal_times <- function(N, times, birth_times, removal_min, removal_max, prob_removal){
  removal_histories <- matrix(0, nrow=N, ncol=length(times))
  for(i in 1:N){
    tmp <- removal_histories[i,]
    
    ## Only find removal time for individual's eligible to be removed 
    if(max(times) >= birth_times[i] + (removal_min)) {
      ## Find removal time
      removal_time <- sample(times[times >= birth_times[i] + (removal_min) & times <= birth_times[i] + (removal_max)], 1)
      tmp[removal_time] <- ifelse(runif(1) < prob_removal, 1, 0) 
      
      #Cannot be removed before birth and removal_min
      tmp[times < birth_times[i] + removal_min] <- NA
      
      #Cannot be removed after removal_max 
      tmp[times > birth_times[i] + removal_max] <- NA
      
      removal_histories[i,] <- tmp
    }
    
    ## For individual's who are never eligible to be removed, add NA for every time step
    if(max(times) <= birth_times[i] + (removal_min)){
      tmp[]<-NA
      removal_histories[i,] <- tmp
      
    }
  }
  
  removal_histories_reshaped <- reshape2::melt(removal_histories) %>% dplyr::mutate(value=as.factor(value))
  colnames(removal_histories_reshaped) <- c("Individual","Time","Removed?")
  removed<-removal_histories_reshaped %>% dplyr::filter(removal_histories_reshaped$`Removed?`==1) %>%  dplyr::mutate(removal_times=Time) %>% dplyr::select(Individual, removal_times) 
  df<- tibble(Individual= 1:N)
  df1<- df %>% dplyr::left_join(removed, by="Individual")
  return(df1$removal_times)
}


#' Generate A Population Demography Data Set 
#' 
#' @description This function generates a tibble of relevant demographic information for each individual in the simulation
#'
#' @param N The number of individuals in the simulation
#' @param times A vector of each time step in the simulation
#' @param birth_times A vector of birth times for each individual; defaults to NULL; if birth_times is not specified then the function will simulate birth times for each individual
#' @param limit A number indicating the youngest age possible by the end of the simulation; defaults to 0 which means individuals can be born up until the second to last time step
#' @param removal_min The minimum age at which an individual can be removed from the population 
#' @param removal_max The maximum age at which an individual can be removed from the population 
#' @param prob_removal The probability that an individual will be removed from the population
#' @param aux A list of the demography columns, the variable options and their distributions; defaults to NULL  
#'
#' @return A tibble of relevant demographic information for each individual in the simulation is returned;  This output matches the required "demography" input for the runserosi function
#' @export
#'
#' @examples generate_pop_demography(10, 1:120, limit=0, removal_min=0, removal_max=120, prob_removal=0.3)
#' 
#' @examples 
#' aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "distribution"=c(0.5,0.5)),"Group"=list("name"="group","options"=c("1", "2", "3", "4"), "distribution"=c(0.25,0.25,0.25,0.25)) )
#' generate_pop_demography(10, 1:120, limit=0, removal_min=0, removal_max=120, prob_removal=0.3, aux=aux)
generate_pop_demography<-function(N, times, birth_times=NULL, limit=0, removal_min, removal_max, prob_removal, aux=NULL){
  if(!is.null(birth_times)){
    birth_tm<-birth_times} 
  else{birth_tm=simulate_birth_times(N, times, limit)}
  
  if(is.null(aux)){
    
    removal_times <- simulate_removal_times(N, times, birth_times=birth_tm, removal_min, removal_max, prob_removal)
    
    
    df<- tibble(
      i=1:N,
      birth= birth_tm,
      removal= removal_times)
    
    exp<- tidyr::expand_grid(1:N, times)
    exp<-dplyr::rename(exp,i="1:N")
    
    dem<- exp %>% dplyr::left_join(df, by="i")
    return(dem)
  }
  
  if(!is.null(aux)){
    removal_times <- simulate_removal_times(N, times, birth_times=birth_tm, removal_min, removal_max, prob_removal)
    
    vars <- NULL
    for(var in seq_along(aux)){
      vars[[var]] <- tibble(
        i=1:N,
        name=unlist(aux[[var]]["name"]),
        value=sample(aux[[var]][["options"]],size=N, prob=aux[[var]][["distribution"]],replace=TRUE)
      )
    }
    vars <- do.call("bind_rows",vars)
    vars <- vars %>% pivot_wider(names_from=name,values_from=value)
    
    df<- tibble(
      i=1:N,
      birth= birth_tm,
      removal= removal_times)
    exp<- tidyr::expand_grid(1:N, times)
    exp<-dplyr::rename(exp,i="1:N")
    dem<- exp %>% dplyr::left_join(df, by="i")
    dem1<- dem %>% left_join(vars, by="i")
    return(dem1)
    
  }
}

#' runserosim function update message
#' 
#' @description This function is used within runserosim to produce messages to inform the user which individual the simulation is currently running
#'
#' @param VERBOSE The interval of individuals in which messages will be printed at
#' @param i The individual that the simulation is currently on 
#'
#' @return A printed statement indicating an individual number is returned 
#' @export
#'
#' @examples
update <- function(VERBOSE, i){ 
  if(!is.null(VERBOSE)){
   if(i %% VERBOSE == 0){
    message(cat("Individual: ", i))
  }
  }
}
