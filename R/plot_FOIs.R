#' Plot The Force Of Infection For All Pathogens
#'
#' @param prob_infections The probability of infection at each time step for all pathogens
#'
#' @return A plot of the force of infection for all pathogens is returned
#' @export
#'
#' @examples
plot_FOIs <- function(prob_infections){
  p_FOI <- ggplot2::ggplot(reshape2::melt(prob_infections) %>% dplyr::mutate(Var2=as.factor(Var2))) + ggplot2::geom_line(ggplot2::aes(x=Var1,y=value,col=Var2)) +
    ggplot2::ylab("Probability of infection") + ggplot2::xlab("Time period") + ggplot2::scale_color_manual(name="Pathogen",values=c("1"="red","2"="orange")) +  ggplot2::theme_bw()
  return(p_FOI)

}


