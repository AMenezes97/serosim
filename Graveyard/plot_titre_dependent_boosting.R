#' Generate A Plot Displaying Proportion Of Boost Received At Each Titre
#'
#' @param start This value is the lower bound of the x axis
#' @param end This value is the upper bound of the x axis
#' @param by This value indicates the increments at which your x axis will be plotted
#' @param kinetics_pars This is your antibody kinetics data frame
#' @param pathogen This is the pathogen you are interested in plotting; 1 or 2
#'
#' @return A plot displaying titre dependent boosting is returned
#' @export
#'
#' @examples
plot_titre_dependent_boosting <- function(start, end, by, kinetics_pars, pathogen){
  if(pathogen==1){
    test_titres <-pmin(seq(start,end,by=by), kinetics_pars$value[18])
    boost_modifier <- (1-kinetics_pars$value[17]*test_titres)
    boost_modifier_dat <- tidyr::tibble(boost_mod=boost_modifier,starting_titres=seq(start,end, by=by))

    g<-ggplot2::ggplot(boost_modifier_dat) + ggplot2::geom_line(ggplot2::aes(x=starting_titres, y=boost_mod)) + ggplot2::theme_bw() +
      ggplot2::xlab("Starting titre") +
      ggplot2::ylab("Proportion of full boost experienced") + ggplot2::scale_y_continuous(limits=c(0,1)) + ggplot2::xlim(start,end)
    return(g)
  }
  if(pathogen==2){
    test_titres <-pmin(seq(start,end,by=by), kinetics_pars$value[36])
    boost_modifier <- (1-kinetics_pars$value[35]*test_titres)
    boost_modifier_dat <- tidyr::tibble(boost_mod=boost_modifier,starting_titres=seq(start,end, by=by))

    g<-ggplot2::ggplot(boost_modifier_dat) + ggplot2::geom_line(ggplot2::aes(x=starting_titres, y=boost_mod)) + ggplot2::theme_bw() +
      ggplot2::xlab("Starting titre") +
      ggplot2::ylab("Proportion of full boost experienced") + ggplot2::scale_y_continuous(limits=c(0,1)) + ggplot2::xlim(start,end)
    return(g)
  }
}


