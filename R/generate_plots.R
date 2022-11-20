#' Plot Titer Mediated Protection Graphs For Each Pathogen
#'
#' @param titer_range The range of possible titers an individual can have at exposure
#' @param titer_prot_midpoint The titer value at which you are 50% protected from infection
#' @param titer_prot_width Determines the chape of the curve
#'
#' @return A plot of the probability of infection given an individual's titer at exposure
#' @export
#'
#' @examples
#' plot_titer_mediated_protection(0:10,5,0.9)
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
  
#' Generate A Plot Displaying Proportion Of Full Boost Received At Each Starting Titer
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
#' plot_titer_dependent_boosting(0,10,1,4,.14)
plot_titer_dependent_boosting <- function(start, end, by, titer_ceiling_threshold, titer_ceiling_gradient){
    test_titres <-pmin(seq(start,end,by=by), titer_ceiling_threshold)
    boost_modifier <- (1-titer_ceiling_gradient*test_titres)
    boost_modifier_dat <- tidyr::tibble(boost_mod=boost_modifier,starting_titers=seq(start,end, by=by))
    
    g<-ggplot2::ggplot(boost_modifier_dat) + ggplot2::geom_line(ggplot2::aes(x=starting_titers, y=boost_mod)) + ggplot2::theme_bw() +
      ggplot2::xlab("Starting titer") +
      ggplot2::ylab("Proportion of full boost experienced") + ggplot2::scale_y_continuous(limits=c(0,1)) + ggplot2::xlim(start,end)
    return(g)
}
  
#' Plot Individual Exposure Probabilities Across Time For All Individuals And Exposure IDs
#'
#' @param exposure_probabilities_long The reshaped data set containing exposure probability for individuals at all time steps for each exposure ID
#'
#' @return A plot of individual exposure probabilities across time for all exposure IDs is returned
#' @export
#'
#' @examples
#' plot_exposure_prob(example_exposure_probabilities)
plot_exposure_prob<-function(exposure_probabilities_long){
  p <- ggplot2::ggplot(exposure_probabilities_long) + 
    ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~x,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual Exposure Probability",
         x="Time",
         y="Individual",
         fill="Probability")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    theme(legend.position="bottom")
  return(p)
}

#' Plot Force of Exposure for Each Individual Across Time For All Individuals And Exposure IDs
#'
#' @param exposure_force_long The reshaped data set containing exposure probability for individuals at all time steps for each exposure ID
#'
#' @return A plot of force of exposure for across time for all individuals and exposure IDs is returned
#' @export
#'
#' @examples
#' plot_exposure_force(example_force_long)
plot_exposure_force<-function(exposure_force_long){
  p <- ggplot2::ggplot(exposure_force_long) + 
    ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~x,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual Force of Exposure",
                  x="Time",
                  y="Individual",
                  fill="Probability")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    theme(legend.position="bottom")
  return(p)
}
#' Plot Individual Exposure Histories
#'
#' @param expsoure_histories The reshaped data set containing exposure history for individuals at all time steps for each exposure ID
#'
#' @return A plot of individual exposures histories across time for all individuals and exposures is returned
#' @export
#'
#' @examples
#' plot_exposure_histories(example_exposure_histories)
plot_exposure_histories <- function(exposure_histories){
  exposure_histories$value <- ifelse(!is.na(exposure_histories$value),
                                     ifelse(exposure_histories$value==1,"Succesful Exposure","No Exposure"), "NA")
  
  p <- ggplot2::ggplot(exposure_histories) + ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + ggplot2::facet_wrap(~x,nrow=2) + ggplot2::theme_bw() + ggplot2::scale_fill_viridis_d() + ggplot2::scale_x_continuous(expand=c(0,0)) + ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual Exposure History",
                  x="Time",
                  y="Individual")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    ggplot2::guides(fill=guide_legend(title="Key")) +
    theme(legend.position="bottom")
  return(p)
}

#' Plot Titers Across Time For All Individuals And biomarkers
#'
#' @param titers The reshaped data set containing antibody titer for individuals at all time steps for each biomarker ID
#'
#' @return A plot of titers across all time steps for all individuals and biomarkers is returned
#' @export
#'
#' @examples
#' plot_titers(example_antibody_states)
plot_titers<- function(titers){
  p <- ggplot2::ggplot(titers) + 
    ggplot2::geom_tile(aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~b,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="True Antibody Titers",
                  x="Time",
                  y="Individual",
                  fill="Titer")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) 
  #+theme(legend.position="bottom")
  
  return(p)
}

#' Plot Observed Titers For One Observation Time 
#' 
#' @description This function should be used when there was only one time step in which titers were observed
#'
#' @param observed_titers The reshaped data set containing observed antibody titers for individuals at all time steps for each biomarker
#'
#' @return A plot of observed titers for all individuals and biomarkers is returned
#' @export
#'
#' @examples
#' plot_obs_titers_one_sample(example_observed_antibody_states)
plot_obs_titers_one_sample<-function(observed_titers){
  p<-ggplot2::ggplot(observed_titers  %>% filter(!is.na(observed))) +
  ggplot2::geom_jitter(ggplot2::aes(x=b, y=observed),
                       height=0,width=0.25) +
  ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
  ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    ggplot2::labs(title="Observed Antibody Titers",
                  x="Biomarker",
                  y="Observed Titer") + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    suppressWarnings(ggplot2::scale_x_discrete(name ="Biomarker", 
                     limits=c(unique(observed_titers$b)), expand = c(0.1, 0.1))) +
    theme(legend.position="bottom")
return(p)
}

#' Plot Observed Titers For Multiple Observation Times and Paired Samples
#' 
#' @description This function should be used when there were multiple time step in which titers were observed
#'
#' @param observed_titers The reshaped data set containing observed antibody titers for individuals at all time steps for each biomarker
#'
#' @return A plot of observed titers for all individuals and biomarkers is returned
#' @export
#'
#' @examples 
plot_obs_titers_paired_sample<-function(observed_titers){
p<- ggplot2::ggplot(observed_titers, aes(x = t, y = observed, group = i)) + 
  ggplot2::geom_line() + 
  ggplot2:: geom_point(size = 2, aes(color=factor(t))) + 
  ggplot2:: facet_wrap(~ b) +
  ggplot2::scale_x_discrete("") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top", 
        panel.grid = element_blank(),
        axis.line.y = element_line(size = .5)) +
  ggplot2::labs(title="Observed Paired Antibody Titers",
                x="Biomarker",
                y="Observed Titer") + 
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
  ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  ggplot2::scale_colour_discrete(name="Sampling Time")
return(p)
}



#' Plot Antibody States and Exposure Histories For A Subset Of Individuals
#'
#' @param titers The reshaped data set containing antibody titer for individuals at all time steps for each biomarker ID
#' @param exposure_histories The reshaped data set containing exposure history for individuals at all time steps for each exposure ID
#' @param subset The number of individuals you want to plot
#' @param demography Tibble of removal time for each individual
#'
#' @return
#' @export
#'
#' @examples
#' plot_subset_individuals_history(example_antibody_states,example_exposure_histories,3,example_demography)
plot_subset_individuals_history <- function(titers, exposure_histories, subset, demography){
  exposure_histories$x <- paste0("Exposure: ", exposure_histories$x)
  titers$b <- paste0("Biomarker: ", titers$b)
  
  exposure_histories_subset<-exposure_histories %>% drop_na() %>% filter(value==1)
  removal_subset <- demography %>% filter(times==1)
  
  sample_indivs <- sample(1:max(demography$i), size=subset)
  
  g<-  ggplot() +
    geom_vline(data=exposure_histories_subset %>% filter(i %in% sample_indivs), aes(xintercept=t, colour=x),linetype="dotted") +
    # geom_vline(data=removal_subset %>% filter(i %in% sample_indivs), aes(xintercept=removal, color="Removal Time"),linetype="solid") +
    geom_line(data=titers %>% filter(i %in% sample_indivs), aes(x=t,y=value,colour=b)) +
    facet_wrap(~i) + theme_bw() +
    scale_color_hue("Key", guide=guide_legend(order=3)) +
    ggplot2::labs(title="Individual Antibody Kinetics",
                  x="Time",
                  y="Antibody Titer") + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
    ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
    ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
    ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
  return(g)
}



