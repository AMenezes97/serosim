#' Plot biomarker mediated protection graphs for each pathogen
#'
#' @param biomarker_range The range of possible biomarker quantities an individual can have at exposure
#' @param biomarker_prot_midpoint The biomarker quantity at which you are 50% protected from infection
#' @param biomarker_prot_width Determines the shape of the curve
#'
#' @return A plot of the probability of infection given an individual's biomarker quantity at exposure
#' @importFrom dplyr tibble
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @export
#'
#' @examples
#' plot_biomarker_mediated_protection(0:10,5,0.9)
plot_biomarker_mediated_protection <- function(biomarker_range, biomarker_prot_midpoint, biomarker_prot_width){
 
  p_infection <- function(phi, biomarker_quantity, alpha1, beta1){
    p <- phi*(1-biomarker_protection(biomarker_quantity, alpha1 , beta1))
    p
  }
  
    #create a data frame with the probability of infection at each biomarker quantity 
    prob_infection <- tibble(biomarker=biomarker_range, prob_infection=p_infection(1, biomarker_range,biomarker_prot_midpoint,biomarker_prot_width))
    #plot probability of infection given biomarker quantity at exposure
    p1<- ggplot2::ggplot(prob_infection) + ggplot2::geom_line(ggplot2::aes(x=biomarker,y=prob_infection)) + ggplot2::theme_bw() + ggplot2::ylab("Probability of infection \n (relative to biomarker quantity of 0)") + ggplot2::xlab("Biomarker quantity at exposure")
    return(p1)
  }
  
#' Generate a plot displaying proportion of full boost received at each starting biomarker quantity 
#'
#' @param start Lower bound of the x axis
#' @param end Upper bound of the x axis
#' @param by Increments at which your x axis will be plotted
#' @param biomarker_ceiling_threshold The maximum biomarker level an individual can have before their boost is limited in size
#' @param biomarker_ceiling_gradient (1-A)/B; Where A is the proportion of the full boost received at or above the biomarker_ceiling_threshold (B)
#'
#' @return A plot displaying biomarker quantity dependent boosting is returned
#' @importFrom dplyr tibble
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 xlim
#' @export
#'
#' @examples
#' plot_biomarker_dependent_boosting(0,10,1,4,.14)
plot_biomarker_dependent_boosting <- function(start, end, by, biomarker_ceiling_threshold, biomarker_ceiling_gradient){
    test_titres <-pmin(seq(start,end,by=by), biomarker_ceiling_threshold)
    boost_modifier <- (1-biomarker_ceiling_gradient*test_titres)
    boost_modifier_dat <- tidyr::tibble(boost_mod=boost_modifier,starting_biomarkers=seq(start,end, by=by))
    
    g<-ggplot2::ggplot(data=boost_modifier_dat) + ggplot2::geom_line(ggplot2::aes(x=starting_biomarkers, y=boost_mod)) + ggplot2::theme_bw() +
      ggplot2::xlab("Starting biomarker") +
      ggplot2::ylab("Proportion of full boost experienced") + ggplot2::scale_y_continuous(limits=c(0,1)) + ggplot2::xlim(start,end)
    return(g)
}
  

#' Plot probability of a successful exposure event for all individuals and exposure events
#'
#' @param exposure_probabilities_long The reshaped data set containing the probability of a successful exposure event for individuals at all time steps for each exposure event
#'
#' @return A plot of the probability of a successful exposure event across time for all individuals and exposure events is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @export
#'
#' @examples
#' plot_exposure_prob(example_exposure_probabilities)
plot_exposure_prob<-function(exposure_probabilities_long){
  p <- ggplot2::ggplot(data=exposure_probabilities_long) + 
    ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~x,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual exposure probability",
         x="Time",
         y="Individual",
         fill="Probability")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                   axis.text.x = element_text(vjust=0.6, size= 8),
                   axis.text.y = element_text(vjust=0.6, size= 8),
                   axis.title.y = element_text(vjust=0.6, size= 10),
                   axis.title.x = element_text(vjust=0.6, size= 10)) +
    theme(legend.position="bottom")
  return(p)
}

#' Plot individual exposure probabilities across time for all individuals and exposure events
#'
#' @param exposure_force_long The reshaped data set containing the probability of an exposure event for individuals at all time steps for each exposure event
#'
#' @return A plot of individual exposure probabilities across time for all exposure events is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom ggplot2 element_text
#' @export
#'
#' @examples
#' plot_exposure_force(example_exposure_force)
plot_exposure_force<-function(exposure_force_long){
  p <- ggplot2::ggplot(exposure_force_long) + 
    ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~x,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual force of exposure",
                  x="Time",
                  y="Individual",
                  fill="Probability")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                   axis.text.x = element_text(vjust=0.6, size= 8),
                   axis.text.y = element_text(vjust=0.6, size= 8),
                   axis.title.y = element_text(vjust=0.6, size= 10),
                   axis.title.x = element_text(vjust=0.6, size= 10)) +
    theme(legend.position="bottom")
  return(p)
}
#' Plot individual immune histories
#'
#' @param immune_histories The reshaped data set containing immune history for individuals at all time steps for each exposure event
#'
#' @return A plot of individual immune histories across time for all individuals and exposure events is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_fill_viridis_d
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 guide_legend
#' @export
#'
#' @examples
#' plot_immune_histories(example_immune_histories)
plot_immune_histories <- function(immune_histories){
  immune_histories$value <- ifelse(!is.na(immune_histories$value),
                                     ifelse(immune_histories$value==1,"Successful exposure","No Exposure"), "NA")
  
  p <- ggplot2::ggplot(immune_histories) + 
    ggplot2::geom_tile(ggplot2::aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~x,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_d() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="Individual immune history",
                  x="Time",
                  y="Individual")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                   axis.text.x = element_text(vjust=0.6, size= 8),
                   axis.text.y = element_text(vjust=0.6, size= 8),
                   axis.title.y = element_text(vjust=0.6, size= 10),
                   axis.title.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::guides(fill=guide_legend(title="Key")) +
    theme(legend.position="bottom")
  return(p)
}

#' Plot biomarker quantities across time for all individuals and biomarkers
#'
#' @param biomarker_states The reshaped data set containing biomarker quantity for individuals at all time steps for each biomarker 
#'
#' @return A plot of biomarker quantities across all time steps for all individuals and biomarkers is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @export
#'
#' @examples
#' plot_biomarker_quantity(example_biomarker_states)
plot_biomarker_quantity<- function(biomarker_states){
  p <- ggplot2::ggplot(data=biomarker_states) + 
    ggplot2::geom_tile(aes(x=t,y=i,fill=value)) + 
    ggplot2::facet_wrap(~b,nrow=2) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::scale_x_continuous(expand=c(0,0)) + 
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::labs(title="True biomarker Quantity",
                  x="Time",
                  y="Individual",
                  fill="Biomarker quantity")   + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                   axis.text.x = element_text(vjust=0.6, size= 8),
                   axis.text.y = element_text(vjust=0.6, size= 8),
                   axis.title.y = element_text(vjust=0.6, size= 10),
                   axis.title.x = element_text(vjust=0.6, size= 10),
                   legend.position="bottom") 
  #+theme(legend.position="bottom")
  
  return(p)
}

#' Plot Observed Biomarker Quantities For One Observation Time 
#' 
#' @description This function should be used when there was only one time step in which biomarker quantities were observed
#'
#' @param observed_biomarker_states The reshaped data set containing observed biomarker quantities for individuals at all time steps for each biomarker
#' @param add_boxplot if TRUE, adds a boxplot to the jittered biomarker quantities with median and 75% percentiles
#' @param add_density if TRUE, adds a smoothed density (violin) plot to the jittered biomarker quantities with 75% quantiles
#'
#' @return A plot of observed biomarker quantities for all individuals and biomarkers is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_violin
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' plot_obs_biomarkers_one_sample(example_observed_biomarker_states)
plot_obs_biomarkers_one_sample<-function(observed_biomarker_states, add_boxplot=FALSE, add_density=FALSE){
  ## Plot last observation for each individual
  observed_biomarker_states <- observed_biomarker_states %>% group_by(i) %>% filter(t == max(t)) %>% ungroup() 
  p<-ggplot2::ggplot(observed_biomarker_states  %>% filter(!is.na(observed)))
  if(add_boxplot){
    p <- p + geom_boxplot(aes(x=b,y=observed,group=b),fill="grey90")
  }
  if(add_density){
    p <- p + geom_violin(aes(x=b,y=observed,group=b),fill="grey90",draw_quantiles=c(0.25,0.5,0.75))
  }
  
  p <- p +
  ggplot2::geom_jitter(ggplot2::aes(x=b, y=observed),
                       height=0,width=0.25) +
  ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                   axis.text.x = element_text(vjust=0.6, size= 8),
                   axis.text.y = element_text(vjust=0.6, size= 8),
                   axis.title.y = element_text(vjust=0.6, size= 10),
                   axis.title.x = element_text(vjust=0.6, size= 10)) +
    ggplot2::labs(title="Observed biomarker quantities",
                  x="Biomarker",
                  y="Observation") + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    suppressWarnings(ggplot2::scale_x_discrete(name ="Biomarker", 
                     limits=c(unique(observed_biomarker_states$b)), expand = c(0.1, 0.1))) +
    theme(legend.position="bottom")
return(p)
}

#' Plot biomarker quantities for multiple observation times and paired samples
#' 
#' @description This function should be used when there were multiple time step in which biomarker quantities were observed
#'
#' @param observed_biomarker_states The reshaped data set containing observed biomarker quantities for individuals at all time steps for each biomarker
#' @param discretize_times if TRUE, plot sampling times as discrete factors. Otherwise, plots sampling time on the x-axis as its raw value.
#'
#' @return A plot of observed biomarker quantities for all individuals and biomarkers is returned
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_colour_discrete
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom dplyr %>%
#' @export
#'
#' @examples 
#' library(dplyr)
#' example_biomarker_states$observed <- example_biomarker_states$value
#' example_biomarker_states_subset <- example_biomarker_states %>% dplyr::filter(t %in% c(1,120))
#' plot_obs_biomarkers_paired_sample(example_biomarker_states_subset)
plot_obs_biomarkers_paired_sample<-function(observed_biomarker_states, discretize_times=TRUE){
  ## Find last 2 times for each individual
  observed_biomarker_states <- observed_biomarker_states %>% dplyr::left_join(
    observed_biomarker_states %>% 
      dplyr::select(i, t) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(i) %>% 
      dplyr::arrange(-t) %>% 
      dplyr::mutate(n_samp=1:n()),
    by=c("i","t")
    ) %>%
    dplyr::filter(n_samp %in% 1:2) %>% 
    dplyr::mutate(n_samp = dplyr::if_else(n_samp == 2, 1, 2))
 

  observed_biomarker_states$biomarker_id <- factor(paste0("Biomarker: ", observed_biomarker_states$b), levels=paste0("Biomarker: ", unique(observed_biomarker_states$b)))
  
  if(discretize_times){
    observed_biomarker_states$t <- observed_biomarker_states$n_samp
    x_breaks <- seq(from=1, to=2,length.out=2)
  } else {
    time_range <- range(observed_biomarker_states$t,na.rm=TRUE)
    x_breaks <- seq(from=time_range[1], to=time_range[2],length.out=5)
  }

p<- ggplot2::ggplot(observed_biomarker_states, aes(x = t, y = observed, group = i)) + 
  ggplot2::geom_line(linewidth=0.25) + 
  ggplot2:: geom_point(size = 2, aes(color=factor(n_samp))) + 
  ggplot2:: facet_wrap(~ biomarker_id) +
  ggplot2::scale_x_continuous("Sampling time", breaks=x_breaks,labels=x_breaks) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        axis.line.y = element_line(linewidth = .5)) +
  ggplot2::labs(title="Paired biomarker quantities",
                y="Observation") + 
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                 axis.text.x = element_text(vjust=0.6, size= 8),
                 axis.text.y = element_text(vjust=0.6, size= 8),
                 axis.title.y = element_text(vjust=0.6, size= 10),
                 axis.title.x = element_text(vjust=0.6, size= 10)) +
  ggplot2::scale_colour_discrete(name="Sample")
return(p)
}



#' Plot biomarker states and immune histories for a subset of individuals
#'
#' @param biomarker_states The reshaped data set containing biomarker quantities for individuals at all time steps for each biomarker 
#' @param immune_histories The reshaped data set containing immune history for individuals at all time steps for each exposure event
#' @param subset The number of individuals you want to plot
#' @param demography Tibble of removal time for each individual
#' @param removal Set to TRUE if individuals are removed during the simulation and removal time is present in demogrpahy; defaults to FALSE
#' @param heatmap if TRUE, returns a heatmap of all biomarker states over time. Overwise plots each biomarker as a line.
#'
#' @return A plot of  biomarker states and exposure histories for a subset of individuals is returned
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_color_hue
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 scale_color_manual
#' @export
#'
#' @examples
#' plot_subset_individuals_history(example_biomarker_states,example_immune_histories,
#' 3,example_demography)
plot_subset_individuals_history <- function(biomarker_states, immune_histories, subset, demography, removal=FALSE, 
                                            heatmap=FALSE){
  immune_histories$exposure_id <- factor(paste0("Exposure: ",immune_histories$x), levels=paste0("Exposure: ", unique(immune_histories$x)))#paste0("Exposure: ", immune_histories$x)
  biomarker_states$biomarker_id <- factor(paste0("Biomarker: ", biomarker_states$b), levels=paste0("Biomarker: ", unique(biomarker_states$b)))
  
  n_exposure_ids <- length(unique(immune_histories$exposure_id))
  n_biomarker_ids <- length(unique(biomarker_states$biomarker_id))
  
  biomarker_colors <- viridis::magma(n_biomarker_ids)
  names(biomarker_colors) <- unique(biomarker_states$biomarker_id)
  exposure_colors <- viridis::viridis(n_exposure_ids)
  names(exposure_colors) <- unique(immune_histories$exposure_id)
  
  all_colors <- c(biomarker_colors, exposure_colors)
  
  immune_histories_subset<-immune_histories %>% drop_na() %>% filter(value==1)
  sample_indivs <- sample(1:max(demography$i), size=subset)
  
  removal_subset <- demography %>% filter(times==1)
  
  g <- ggplot()
  
  if(removal==TRUE){
    g <- g + geom_vline(data=removal_subset %>% filter(i %in% sample_indivs), 
                        aes(xintercept=removal, color="Removal time"),linetype="solid")
    all_colors <- c(all_colors, "Removal time" ="red")
  }

  if(heatmap){
    g <- g + geom_tile(data=biomarker_states %>% filter(i %in% sample_indivs), aes(x=t,y=biomarker_id,fill=value)) +
      scale_fill_viridis_c(name="Biomarker quantity", option="magma")
    
  } else {
    g <- g +
      geom_line(data=biomarker_states %>% filter(i %in% sample_indivs), aes(x=t,y=value,colour=biomarker_id))
  }
  g <- g + geom_vline(data=immune_histories_subset %>% filter(i %in% sample_indivs), 
                      aes(xintercept=t, colour=exposure_id),linetype="dashed",linewidth=0.75)
  g <- g +
      facet_wrap(~i) + 
      theme_bw() +
      #scale_color_hue("Biomarker and Exposure Key", guide=guide_legend(order=3)) +
      #ggplot2::scale_color_viridis_d(name="Biomarker and Exposure Key") + 
    ggplot2::scale_color_manual(name="Biomarker and exposure key",values=all_colors)+
    guides(color=guide_legend(ncol=4)) +
      ggplot2::labs(title="Individual biomarker kinetics",
                    x="Time",
                    y="Biomarker quantity") + 
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                     axis.text.x = element_text(vjust=0.6, size= 8),
                     axis.text.y = element_text(vjust=0.6, size= 8),
                     axis.title.y = element_text(vjust=0.6, size= 10),
                     axis.title.x = element_text(vjust=0.6, size= 10)) +
      theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
  
    
  return(g)
}

#' Plots example trajectories of the provided antibody kinetics model
#' 
#' @param N Number of trajectories to simulate
#' @param times Vector of times to solve model over
#' @param draw_parameters_fn Pointer to function used for drawing random kinetics parameters, see \code{\link{draw_parameters_fixed_fx}}
#' @inheritParams runserosim
#' @return A ggplot2 object
#' @importFrom dplyr filter
#' @importFrom dplyr group_by 
#' @importFrom dplyr summarize 
#' @importFrom dplyr pull
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_fill_viridis_d
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom stats complete.cases
#' @export
#' @examples
#' library(dplyr)
#' model_pars <- reformat_biomarker_map(example_model_pars_biphasic)%>% tidyr::drop_na()
#' draw_parameters_random_fx(1,1,1,NULL,NULL,model_pars)
#' biomarker_map <- dplyr::tibble(exposure_id=c(1,1,2),biomarker_id=c(1,2,1))
#' plot_antibody_model(antibody_model_biphasic, 50, model_pars=model_pars,
#' draw_parameters_fn = draw_parameters_random_fx, biomarker_map=biomarker_map)
#' plot_antibody_model(antibody_model_biphasic, 50, model_pars=model_pars,
#' draw_parameters_fn = draw_parameters_fixed_fx, biomarker_map=biomarker_map)
plot_antibody_model <- function(antibody_model,N=100, times=seq(1,50,by=1),model_pars,biomarker_map, 
                                demography=NULL, draw_parameters_fn=draw_parameters_fixed_fx, ...){
  exposure_ids <- unique(biomarker_map$exposure_id)
  biomarker_ids <- unique(biomarker_map$biomarker_id)
  
    ## Go through for all times and plot random trajectories
    indivs <- 1:N
    immune_histories_tmp <- array(0, dim=c(N, length(times), length(exposure_ids)))
    antibody_states_all <- list()
    for(x in exposure_ids){
        biomarkers_tmp <- biomarker_map %>% filter(exposure_id == x) %>% pull(biomarker_id)
        n_biomarkers_tmp <- length(biomarkers_tmp)
        antibody_states <- array(0, dim=c(N, length(times),length(biomarker_ids)))
        kinetics_pars_tmp <- list()
        for(i in indivs){
            immune_histories_tmp[i,1,x] <- 1
            for(b in biomarker_ids){
                if(b %in% biomarkers_tmp){
                    kinetics_pars_tmp <- list(draw_parameters_fn(i, 1, x, demography,antibody_states, model_pars, ...))
                    kinetics_pars_tmp[[1]] <- kinetics_pars_tmp[[1]][complete.cases(kinetics_pars_tmp[[1]]),]
                    antibody_states[i,,b] <- vapply(times,function(t) antibody_model(1, t, b, immune_histories_tmp,antibody_states, 
                                                                                     kinetics_pars_tmp, biomarker_map, ...), numeric(1))
                } else {
                    antibody_states[i,,b] <- NA
                }
            }
        }
        antibody_states <- reshape2::melt(antibody_states) %>% drop_na()
        colnames(antibody_states) <- c("i","t","b","titer")
        antibody_states$x <- x
        antibody_states_all[[x]] <- antibody_states
    }
    antibody_states_all <- do.call("bind_rows", antibody_states_all)
    antibody_states_all$exposure_id <-factor(paste0("Exposure: ",antibody_states_all$x), levels=paste0("Exposure: ", unique(antibody_states_all$x)))
    antibody_states_all$biomarker_id <- factor(paste0("Biomarker: ", antibody_states_all$b), levels=paste0("Biomarker: ", unique(antibody_states_all$b)))
    
    antibody_states_summ <- antibody_states_all %>% group_by(t,b,x) %>% summarize(mean_titer=mean(titer))
    antibody_states_summ$biomarker_id <- factor(paste0("Biomarker: ", antibody_states_summ$b), levels=paste0("Biomarker: ", unique(antibody_states_summ$b)))
    antibody_states_summ$exposure_id <-factor(paste0("Exposure: ",antibody_states_summ$x), levels=paste0("Exposure: ", unique(antibody_states_summ$x)))    
    ## Create biomarker colors
    n_biomarker_ids <- length(unique(antibody_states_all$biomarker_id))
    biomarker_colors <- viridis::magma(n_biomarker_ids)
    names(biomarker_colors) <- unique(antibody_states_all$biomarker_id)

    p <-  ggplot2::ggplot(antibody_states_all) +
      ggplot2::geom_line(aes(x=t,y=titer,col=biomarker_id,group=i),alpha=0.25) +
      ggplot2::geom_line(data=antibody_states_summ,aes(x=t,y=mean_titer,col=biomarker_id),linewidth=1) +
        ggplot2::theme_bw() +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                     axis.text.x = element_text(vjust=0.6, size= 8),
                     axis.text.y = element_text(vjust=0.6, size= 8),
                     axis.title.y = element_text(vjust=0.6, size= 10),
                     axis.title.x = element_text(vjust=0.6, size= 10)) +
      ggplot2::xlab("Time since infection") +
      ggplot2::ylab("Biomarker quantity") +
        ggplot2::scale_color_manual(name="Biomarker",values=biomarker_colors) +
      ggplot2::facet_grid(b~exposure_id,scales="free_y")
    return(p)
    
}

#' Plots the probability of exposure over time for the provided exposure models
#' 
#' @inheritParams runserosim
#' @param indivs (optional) vector of individuals to plot exposure probabilities for. This is important if the `demography` table contains information on more individuals than you wish to plot
#' @param times Vector of times to solve model over
#' @param n_groups Number of groups corresponding to \code{foe_pars}
#' @param n_exposures Number of exposure types corresponding to \code{foe_pars}
#' @param foe_pars Generic object containing all parameters needed to solve \code{exposure_model}
#' @return A ggplot2 object
#' @importFrom dplyr filter
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom dplyr %>% 
#' @export
#' @examples 
#' ## Basic exposure model with demography modifier
#' times <- seq(1,120,1)
#' n_groups <- 1
#' n_exposures <- 2
#' foe_pars <- array(NA, dim=c(n_groups,length(times),n_exposures))
#' foe_pars[1,,1] <- 0.01
#' foe_pars[1,,2] <- 0.005 
#' aux <- list("SES"=list("name"="SES","options"=c("low","high"), "distribution"=c(0.5,0.5)))
#' demography <- generate_pop_demography(N=5, times, age_min=0, removal_min=1, 
#' removal_max=120, prob_removal=0.2, aux=aux)
#' dem_mod <- dplyr::tibble(exposure_id=c(1,1,2,2),column=c("SES","SES","SES","SES"),
#'                  value=c("low","high","low","high"),modifier=c(1,0.75,1,0.5))
#' 
#' plot_exposure_model(indivs=1:5, exposure_model=exposure_model_dem_mod, 
#' times=times,1,2,foe_pars=foe_pars,demography = demography,dem_mod=dem_mod)
#'                     
#' ## SIR model with two groups and two exposure types                 
#' foe_pars <- dplyr::bind_rows(
#'                       dplyr::tibble(x=1,g=1,name=c("beta","gamma","I0","R0","t0"),
#'                       value=c(0.3,0.2,0.00001,0,0)),
#'                       dplyr::tibble(x=2,g=1,name=c("beta","gamma","I0","R0","t0"),
#'                       value=c(0.35,0.25,0.00001,0,200)),
#'                       dplyr::tibble(x=1,g=2,name=c("beta","gamma","I0","R0","t0"),
#'                       value=c(0.5,0.2,0.00005,0,0)),
#'                       dplyr::tibble(x=2,g=2,name=c("beta","gamma","I0","R0","t0"),
#'                       value=c(0.27,0.2,0.00001,0,50))
#'                       )
#' plot_exposure_model(exposure_model=exposure_model_sir, times=seq(1,365,by=1),
#' n_groups = 2,n_exposures = 2,foe_pars=foe_pars)
#'                     
plot_exposure_model <- function(indivs=1, exposure_model, times, n_groups=1, n_exposures=1, foe_pars, demography=NULL, ...){
    n_times <- length(times)
    
    ## Solve exposure probability for each demographic element
    if(!is.null(demography)){
        demography <- demography %>% filter(i %in% indivs)
        n_indivs <- length(unique(demography$i))
    } else {
        n_indivs <- 1
    }
    ## Solve exposure probability for each exposure type, for each group, for each time
    foe_all <- NULL
    for(i in 1:n_indivs){
        foe <- array(NA, dim=c(n_groups,n_times,n_exposures))
        for(g in 1:n_groups){
            for(x in 1:n_exposures){
                foe[g,,x] <- unlist(vapply(times, function(t) exposure_model(i, t, x, g, foe_pars, demography, ...),numeric(1)))
            }
        }
        foe <- reshape2::melt(foe)
        colnames(foe) <- c("Group","Time","Exposure ID","value")
        foe$`Exposure ID` <- paste0("Exposure ID: ", foe$`Exposure ID`)
        foe$Individual <- paste0("Individual: ",i)
        foe_all[[i]] <- foe
    }
    foe_all <- do.call("bind_rows", foe_all)
    foe_all$Group <- as.factor(foe_all$Group)
    foe_all$`Exposure ID` <- as.factor(foe_all$`Exposure ID`)
    p <- ggplot(data=foe_all) + 
        geom_line(aes(x=Time,col=Group,y=value,group=interaction(Group,`Exposure ID`))) + 
        ylab("Probability of exposure per unit time") +
        xlab("Time") +
        theme_bw() +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, size=12),
                     axis.text.x = element_text(vjust=0.6, size= 8),
                     axis.text.y = element_text(vjust=0.6, size= 8),
                     axis.title.y = element_text(vjust=0.6, size= 10),
                     axis.title.x = element_text(vjust=0.6, size= 10)) +
        facet_grid(Individual~`Exposure ID`)
    return(p)
}
