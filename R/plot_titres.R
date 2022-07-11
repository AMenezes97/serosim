#' Plot Titres Across Time For All Individuals And Pathogens
#'
#' @param titres The reshaped data set containing antibody titre for individuals at all time steps for each pathogen
#'
#' @return A plot of titres across time for all individuals and pathogens is returned
#' @export
#'
#' @examples
plot_titres<- function(titres){
  p <- ggplot(titres) + geom_tile(aes(x=Time,y=Individual,fill=Titre)) + facet_wrap(~Pathogen,nrow=2) + theme_bw() + scale_fill_viridis_c() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

    return(p)
}
