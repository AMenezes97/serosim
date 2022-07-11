#' Plot Infection Probabilities Across Time For All Individuals And Pathogens
#'
#' @param all_infection_probs The reshaped data set containing infection probability for individuals at all time steps for each pathogen
#'
#' @return A plot of infection probabilities across time for all individuals and pathogens is returned
#' @export
#'
#' @examples
plot_all_infection_probs<-function(all_infection_probs){
  p <- ggplot2::ggplot(all_infection_probs) + ggplot2::geom_tile(ggplot2::aes(x=Time,y=Individual,fill=`Infection probability`)) + ggplot2::facet_wrap(~Pathogen,nrow=2) + ggplot2::theme_bw() + ggplot2::scale_fill_viridis_c() + ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))
  return(p)
}
