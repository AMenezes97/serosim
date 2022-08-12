#' Convert Normal Distribution Standard Deviation to Log-Normal Distribution Standard Deviation
#'
#' @param normmean The normal distribution mean
#' @param normsd The normal distribution standard deviation
#'
#' @return The log-normal distribution standard deviation is printed
#' @export
#'
#' @examples normal_to_lognormal_sd(5000,400)
normal_to_lognormal_sd <- function(normmean, normsd) {
  phi <- sqrt(normsd ^ 2 + normmean ^ 2)
  sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
  return(sdlog)
}
