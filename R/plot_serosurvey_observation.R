#' Plot Serosurvey Observation
#'
#' @param titres The data set with individual's titres at each time step
#' @param obs_time The observation time when titres were measured
#' @param titre The titre type to plot, "true" for true titre or "observed" for observed titre
#' @param Events If true, point color indicates true vaccination and infection status; titres data set must be in the full_truth form (see full_truth function); defaults to "false"
#'
#' @return A plot of titres for all individual's at the observation time for all pathogens is returned
#' @export
#'
#' @examples
plot_serosurvey_observation<-function(titres, obs_time, titre="observed", Events="false"){
  if(titre=="true" & Events== "true"){
    p<-ggplot2::ggplot(titres  %>% filter(Time == obs_time)) +
      ggplot2::geom_jitter(ggplot2::aes(x=Pathogen, y=`Titre`,col=Events),
                  height=0,width=0.25)+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 0)) +
      ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
      ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 15))
    return(p)
  }
  if(titre=="observed" & Events== "true"){
    p<-ggplot2::ggplot(titres  %>% filter(Time == obs_time)) +
      ggplot2::geom_jitter(ggplot2::aes(x=Pathogen,y=`Observed titre`,col=Events),
                           height=0,width=0.25)+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 0)) +
      ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
      ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 15))
    return(p)
  }
  if(titre=="true" & Events== "false"){
    p<-ggplot2::ggplot(titres %>% filter(Time == obs_time)) +
      ggplot2::geom_jitter(ggplot2::aes(x=Pathogen,y=`Titre`,col=Pathogen),
                  height=0,width=0.25)+
      ggplot2::scale_color_manual(values=c("1"="red","2"="orange")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 0)) +
      ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
      ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 15))
    return(p)
  }
  if(titre=="observed" & Events =="false"){
    p<- ggplot2::ggplot(titres %>% filter(Time == obs_time)) +
      ggplot2::geom_jitter(ggplot2::aes(x=Pathogen,y=`Observed titre`,col=Pathogen),
                  height=0,width=0.25)+
      ggplot2::scale_color_manual(values=c("1"="red","2"="orange")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 0)) +
      ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
      ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 15))
    return(p)
  }
}
