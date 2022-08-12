#' Generate A Plot Of The Age Distribution
#'
#' @param birth_times A vector of all of the birth times of the individuals
#'
#' @return A histogram plot with the number of individuals born at each time step is returned
#' @export
#'
#' @examples
plot_age_distribution <-function(birth_times){
  g<-ggplot2::ggplot(tidyr::tibble(`Birth times`=birth_times)) +
    ggplot2::geom_histogram(ggplot2::aes(x=`Birth times`), binwidth=1) +
                              ggplot2::theme_bw() + ggplot2::xlab("Birth times (month)") + ggplot2::ylab("Number of individuals")
  return(g)
}
