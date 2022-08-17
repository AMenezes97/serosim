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
  
#' Generate A Plot Displaying Proportion Of FUll Boost Received At Each Starting Titer
#'
#' @param start Lower bound of the x axis
#' @param end Upper bound of the x axis
#' @param by Increments at which your x axis will be plotted
#' @param titer_ceiling_threshold The maximum titer level an individual can have before their boost is limited in size
#' @param titer_ceiling_gradient (1-A)/B; Where A is the proportion of the full boost received at or above the titer_ceiling_threshold (B)
#'
#' @return A plot displaying titer dependent boosting is returned
#' @export
#'
#' @examples
plot_titer_dependent_boosting <- function(start, end, by, titer_ceiling_threshold, titer_ceiling_gradient){
    test_titres <-pmin(seq(start,end,by=by), titer_ceiling_threshold)
    boost_modifier <- (1-titer_ceiling_gradient*test_titres)
    boost_modifier_dat <- tidyr::tibble(boost_mod=boost_modifier,starting_titers=seq(start,end, by=by))
    
    g<-ggplot2::ggplot(boost_modifier_dat) + ggplot2::geom_line(ggplot2::aes(x=starting_titers, y=boost_mod)) + ggplot2::theme_bw() +
      ggplot2::xlab("Starting titer") +
      ggplot2::ylab("Proportion of full boost experienced") + ggplot2::scale_y_continuous(limits=c(0,1)) + ggplot2::xlim(start,end)
    return(g)
}
  
  
  
  
  
  
  
  
  
  








#' Title
#'
#' @param infection_histories The reshaped data set containing infection histories for individuals at all time steps for each pathogen
#'
#' @return A plot of infection probabilities across time for all individuals and pathogens is returned
#' @export
#'
#' @examples
plot_infection_histories <- function(infection_histories){
  p<- ggplot2::ggplot(infection_histories) + ggplot2::geom_tile(ggplot2::aes(x=t,y=i)) + ggplot2::facet_wrap(~e,nrow=2) + ggplot2::theme_bw() + ggplot2::scale_fill_viridis_d() + ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))
  return(p)
}

