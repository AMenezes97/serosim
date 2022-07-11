#' Title
#'
#' @param infection_histories The reshaped data set containing infection histories for individuals at all time steps for each pathogen
#'
#' @return A plot of infection probabilities across time for all individuals and pathogens is returned
#' @export
#'
#' @examples
plot_infection_histories <- function(infection_histories){
  p<- ggplot2::ggplot(infection_histories) + ggplot2::geom_tile(ggplot2::aes(x=Time,y=Individual,fill=`Infected?`)) + ggplot2::facet_wrap(~Pathogen,nrow=2) + ggplot2::theme_bw() + ggplot2::scale_fill_viridis_d() + ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0))
  return(p)
}