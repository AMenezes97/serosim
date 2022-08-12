#' Plot Infection and Vaccine History For A Subset Of Individuals
#'
#' @param infection_histories Infection history data set
#' @param vaccine_histories_reshaped Vaccine history data set
#' @param titres  Titre data set
#' @param subset The number of individuals you want to plot
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return A plot of infection and vaccine history for a subset of individuals is returned
#' @export
#'
#' @examples
plot_subset_individuals_history <- function(infection_histories, vaccine_histories_reshaped, titres, subset, N_pathogens){
  if(N_pathogens==1){
    infection_histories_subset<- infection_histories %>% drop_na()%>% filter(`Infected?`==1)
    vaccine_histories_subset<- vaccine_histories_reshaped %>% drop_na()%>% filter(`Vaccinated?`==1)

    sample_indivs <- sample(1:N, size=subset)

   g<-  ggplot2::ggplot() +
      ggplot2::geom_vline(data=vaccine_histories_subset %>% filter(Individual %in% sample_indivs), ggplot2::aes(xintercept=Time),col="blue",linetype="dotted") +
      ggplot2::geom_vline(data=infection_histories_subset %>% filter(Individual %in% sample_indivs), ggplot2::aes(xintercept=Time,col=Pathogen),linetype="dashed") +
      ggplot2::geom_line(data=titres %>% filter(Individual %in% sample_indivs),ggplot2::aes(x=Time,y=Titre,col=Pathogen)) +
      ggplot2::geom_point(data=titres %>% filter(Individual %in% sample_indivs, Time == obs_time),ggplot2::aes(x=Time,y=`Observed titre`,col=Pathogen),shape=4) +
      ggplot2::scale_color_manual(values=c("1"="red")) +
      ggplot2::facet_wrap(~Individual)
    return(g)
  }
  if(N_pathogens==2){

    infection_histories_subset<- infection_histories %>% drop_na()%>% filter(`Infected?`==1)
    vaccine_histories_subset<- vaccine_histories_reshaped %>% drop_na()%>% filter(`Vaccinated?`==1)

    sample_indivs <- sample(1:N, size=subset)

    g<- ggplot2::ggplot() +
      ggplot2::geom_vline(data=vaccine_histories_subset %>% filter(Individual %in% sample_indivs), ggplot2::aes(xintercept=Time),col="blue",linetype="dotted") +
      ggplot2::geom_vline(data=infection_histories_subset %>% filter(Individual %in% sample_indivs), ggplot2::aes(xintercept=Time,col=Pathogen),linetype="dashed") +
      ggplot2::geom_line(data=titres %>% filter(Individual %in% sample_indivs),ggplot2::aes(x=Time,y=Titre,col=Pathogen)) +
      ggplot2::geom_point(data=titres %>% filter(Individual %in% sample_indivs, Time == obs_time),ggplot2::aes(x=Time,y=`Observed titre`,col=Pathogen),shape=4) +
      ggplot2::scale_color_manual(values=c("1"="red","2"="orange")) +
      ggplot2::facet_wrap(~Individual)
    return(g)
  }
}
