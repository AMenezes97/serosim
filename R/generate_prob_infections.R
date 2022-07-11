#' Generate The Probability Of Infection At Each Time Step For All Pathogens
#'
#' @param FOIs The Force of Infection at each time step for each pathogen; only needed for type=1
#' @param pathogens The number of pathogens
#' @param overall_infection_probs The cumulative probability of infection for each pathogen
#' @param type type=1 if you are simulating random FOIs, type=2 if you are inputting constant force of infection; default=1
#' @param times The number of time steps in the simulation
#' @param set_prob1 The probability of infection for pathogen 1 at all time steps; only needs to be defined for type=2
#' @param set_prob2 The probability of infection for pathogen 2 at all time steps; only needs to be defined for type=2
#'
#' @return The probability of infection at each time step for all pathogens is returned
#' @export
#'
#' @examples
generate_prob_infections <- function(FOIs, pathogens, overall_infection_probs, type=1, times, set_prob1, set_prob2){
  if(type==1){
    #Generates probability of infection given the randomly generated FOI
    prob_infections <- apply(FOIs, pathogens, function(x){
      tmp <- 1/(1+exp(-x))
      tmp <- tmp/sum(tmp)
    })
    prob_infections <- t(t(prob_infections)*(overall_infection_probs=overall_infection_probs))
    return(prob_infections)
  }
  if(type==2){
    #Constant circulation of pathogen 1 and no circulation of pathogen 2
    prob_infections<- matrix(data=rep(0), nrow = max(times=times), ncol=pathogens)
    prob_infections[,1]<-rep(set_prob1) #input probability of infection for pathogen 1
    prob_infections[,2]<-rep(set_prob2) #input probability of infection for pathogen 2
    prob_infections <- t(t(prob_infections)*(overall_infection_probs=overall_infection_probs))
    return(prob_infections)
  }
}
