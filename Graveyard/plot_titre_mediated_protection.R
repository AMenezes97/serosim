#' Plot Titre Mediated Protection Graphs For Each Pathogen
#'
#' @param titre_range The range of possible titres an individual can have at exposure
#' @param N_pathogens The number of pathogens in the simulation
#' @param titre_prot_midpoint The titre value at which you are 50% protected from infection
#' @param titre_prot_width
#'
#' @return A plot of the probability of infection given an individual's titre at exposure for all pathogens is returned
#' @export
#'
#' @examples
plot_titre_mediated_protection <- function(titre_range, N_pathogens, titre_prot_midpoint, titre_prot_width){
  if(N_pathogens==1){
    #create a data frame with the probability of infection at each titre level
    prob_infection <- tibble(titre=titre_range, prob_infection=p_infection(1, titre_range,titre_prot_midpoint,titre_prot_width))
    #plot probability of infection given titre level at exposure
    p1<- ggplot2::ggplot(prob_infection) + ggplot2::geom_line(ggplot2::aes(x=titre,y=prob_infection)) + ggplot2::theme_bw() + ggplot2::ylab("Probability of infection (relative to 0 titre)") + ggplot2::xlab("Titre at exposure")
    return(p1)
  }

  if(N_pathogens==2){

    #create a data frame with the probability of infection at each titre level 1-10000
    prob_infection1 <- tibble(titre=titre_range,prob_infection=p_infection(1, titre_range,titre_prot_midpoint[1],titre_prot_width[1]))
    #plot probability of infection given titre level at exposure for pathogen 1
    p1<- ggplot2::ggplot(prob_infection1) + ggplot2::geom_line(ggplot2::aes(x=titre,y=prob_infection)) + ggplot2::theme_bw() + ggplot2::ylab("Probability of infection to Pathogen 1 (relative to 0 titre)") + ggplot2::xlab("Titre at exposure")

    #create a data frame with the probability of infection at each titre level 1-10000
    prob_infection2 <- tibble(titre=titre_range,prob_infection=p_infection(1, titre_range,titre_prot_midpoint[2],titre_prot_width[2]))
    #plot probability of infection given titre level at exposure
    p2<- ggplot2::ggplot(prob_infection2) + ggplot2::geom_line(ggplot2::aes(x=titre,y=prob_infection)) + ggplot2::theme_bw() + ggplot2::ylab("Probability of infection to Pathogen 2 (relative to 0 titre)") + ggplot2::xlab("Titre at exposure")

    return(p1 + p2 + patchwork::plot_layout(nrow=2,heights=c(5,5)))
  }
}
