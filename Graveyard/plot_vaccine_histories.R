#' Plot Vaccine Histories
#'
#' @param vaccine_histories  A matrix with each individual's vaccination histories
#'
#' @return A plot displaying individual's vaccine histories is returned
#' @export
#'
#' @examples
plot_vaccine_histories <- function(vaccine_histories){
  vaccine_histories_melted <- reshape2::melt(vaccine_histories) %>% mutate(value=as.factor(value))
  #Label the columns
  colnames(vaccine_histories_melted) <- c("Individual","Time","Vaccinated?")
  #Look at vaccine histories as trajectories
  g<-ggplot2::ggplot(vaccine_histories_melted) + ggplot2::geom_tile(ggplot2::aes(x=Time,y=Individual,fill=`Vaccinated?`)) + ggplot2::theme_bw() + xlab("Time") + ylab("Individual")
  return(g)
}
