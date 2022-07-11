#' Plot Serosurvey Titre Distribution
#'
#' @param titresobs The titre data set for all individuals at the observation time
#' @param titre Indicate whether you want to plot "true" or "observed" titres; defaults to "observed"
#' @param seronegative Indicate whether you want to keep or remove seronegative individuals; defaults to "keep"
#' @param bin The size of the bins for the histogram
#' @param color The variable to label the points either by "events" or "pathogen"; defaults to "pathogen"
#'
#' @return A histogram is returned
#' @export
#'
#' @examples
plot_serosurvey_distribution<- function(titresobs, titre="observed", seronegative="keep", bin, color="pathogen"){
  if(titre=="observed" & seronegative=="keep" & color=="pathogen"){
    g<- ggplot(titresobs, aes(x='Observed titre')) +
      geom_histogram(aes(fill=factor(Pathogen)), binwidth =bin) +
      labs(title="Distribution of true titres",
           subtitle="",
           x="Titre",
           y="Number of individuals",
           fill="Pathogen")   +
      theme(plot.title = element_text(hjust = 0.5))  +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(g)
  }

  if(titre=="true" & seronegative=="keep" & color=="pathogen"){
  g <- ggplot2::ggplot(titresobs, ggplot2::aes(Titre))
  p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Pathogen)), binwidth = bin) +
    ggplot2::labs(title="Distribution of observed titres",
         subtitle="",
         x="Titre",
         y="Number of individuals",
         fill="Pathogen")   +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
    ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(p)
  }
  if(titre=="observed" & seronegative=="keep" & color=="events"){
    g <- ggplot2::ggplot(titresobs, ggplot2::aes('Observed titre'))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Events)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Events")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }

  if(titre=="true" & seronegative=="keep" & color=="events"){
    g <- ggplot2::ggplot(titresobs, ggplot2::aes(Titre))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Events)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Events")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }
  if(titre=="observed" & seronegative=="remove" & color=="pathogen"){
    g <- ggplot2::ggplot(titresobs %>%  dplyr::filter('Observed titre'>0), ggplot2::aes('Observed titre'))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Pathogen)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Pathogen")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }

  if(titre=="true" & seronegative=="remove" & color=="pathogen"){
    g <- ggplot2::ggplot(titresobs %>%  dplyr::filter(Titre>0), ggplot2::aes(Titre))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Pathogen)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Pathogen")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }
  if(titre=="observed" & seronegative=="remove" & color=="events"){
    g <- ggplot2::ggplot(titresobs %>%  dplyr::filter('Observed titre'>0), ggplot2::aes('Observed titre'))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Events)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Events")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }

  if(titre=="true" & seronegative=="remove" & color=="events"){
    g <- ggplot2::ggplot(titresobs %>%  dplyr::filter(Titre>0), ggplot2::aes(Titre))
    p <- g + ggplot2::geom_histogram(ggplot2::aes(fill=factor(Events)), binwidth = bin) +
      ggplot2::labs(title="Distribution of observed titres",
                    subtitle="",
                    x="Titre",
                    y="Number of individuals",
                    fill="Events")   +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))  +
      ggplot2:: theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    return(p)
  }

}

