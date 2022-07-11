#' Convert Normal Distribution Mean to Log-Normal Distribution Mean
#'
#' @param normmean The normal distribution mean
#' @param normsd The normal distribution standard deviation
#'
#' @return The log-normal distribution mean is printed
#' @export
#'
#' @examples normal_to_lognormal_mean(5000,400)
normal_to_lognormal_mean <- function(normmean, normsd) {
  phi <- sqrt(normsd ^ 2 + normmean ^ 2)
  meanlog <- log(normmean ^ 2 / phi)
  return(meanlog)
}
