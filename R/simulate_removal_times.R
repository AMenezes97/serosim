#' Simulate  Removal Times for Individuals 
#'
#' @param N The number of individuals in the simulation
#' @param times The total number of time steps in the simulation
#' @param birth_times A vector of all individual's birth times
#' @param removal_min The minimum age at which an individual can be removed from the population 
#' @param removal_max The maximum age at which an individual can be removed from the population 
#' @param prob_removal The probability that an individual will be removed from the population
#'
#' @return A vector of all individual's removal times is returned 
#' @export
#'
#' @examples
simulate_removal_times <- function(N, times, birth_times, removal_min, removal_max, prob_removal){
  removal_histories <- matrix(0, nrow=N, ncol=length(times))
  for(i in 1:N){
    tmp <- removal_histories[i,]
    
    #Cannot be removed before birth and removal_min
    tmp[times < birth_times[i] + removal_min] <- NA
    
    #Find removal time
    if(birth_times[i]>removal_min){
      removal_time <- sample(times[times > birth_times[i] + (removal_min) & times < birth_times[i] + (removal_max)], 1)
      tmp[removal_time] <- ifelse(runif(1) < prob_removal, 1, 0) #this step can be removed?
    }
    
    removal_histories[i,] <- tmp
  }
 
  removal_histories_reshaped <- reshape2::melt(removal_histories) %>% dplyr::mutate(value=as.factor(value))
  colnames(removal_histories_reshaped) <- c("Individual","Time","Removed?")
  removed<-removal_histories_reshaped %>% dplyr::filter(removal_histories_reshaped$`Removed?`==1) %>%  dplyr::mutate(removal_times=Time) %>% dplyr::select(Individual, removal_times) 
  df<- tibble(Individual= 1:N)
  df1<- df %>% dplyr::left_join(removed, by="Individual")
  return(df1$removal_times)
}

