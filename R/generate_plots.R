#' Plot Titer Mediated Protection Graphs For Each Pathogen
#'
#' @param titer_range The range of possible titers an individual can have at exposure
#' @param titer_prot_midpoint The titer value at which you are 50% protected from infection
#' @param titer_prot_width
#'
#' @return A plot of the probability of infection given an individual's titer at exposure
#' @export
#'
#' @examples
plot_titer_mediated_protection <- function(titer_range, titer_prot_midpoint, titer_prot_width){
  ## Create a function to calculate the risk of infection at a given titer
  titer_protection <- function(titer, alpha1, beta1){
    risk <- 1 - 1/(1 + exp(beta1*(titer - alpha1)))
    return(risk)
  }
  
  p_infection <- function(phi, titer, alpha1, beta1){
    p <- phi*(1-titer_protection(titer, alpha1 , beta1))
    p
  }
  
    #create a data frame with the probability of infection at each titre level
    prob_infection <- tibble(titer=titer_range, prob_infection=p_infection(1, titer_range,titer_prot_midpoint,titer_prot_width))
    #plot probability of infection given titre level at exposure
    p1<- ggplot2::ggplot(prob_infection) + ggplot2::geom_line(ggplot2::aes(x=titer,y=prob_infection)) + ggplot2::theme_bw() + ggplot2::ylab("Probability of infection (relative to 0 titre)") + ggplot2::xlab("Titer at exposure")
    return(p1)
  }
  
