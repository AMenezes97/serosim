#' Draw A Random Parameter From A Distribution
#'
#' @param par1 The mean of the distribution
#' @param par2 The standard deviation of the distribution
#' @param distribution The distribution to use; defaults to "log-normal"
#'
#' @return A parameter value gets printed
#' @export
#'
#' @examples draw_parameters(5000,400,"log-normal")
#' @examples draw_parameters(5000,400,"normal")
draw_parameters <- function(par1, par2, distribution="log-normal"){
  if(par1==0 & par2==0){
    return(0)
  }
  if(distribution=="log-normal"){
    mean<-normal_to_lognormal_mean(par1,par2)
    sd<-normal_to_lognormal_sd(par1,par2)
    return(stats::rlnorm(1, mean, sd))
  }
  if(distribution=="normal"){
    return(stats::rnorm(1,par1,par2))
  }
}
