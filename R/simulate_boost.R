#' Simulate Titre Boost
#'
#' @param current_titre An individual's current titre at the time of this new boost
#' @param boost_mean The mean of the distribution which the boost is selected from
#' @param boost_sd The standard deviation of the distribution which the boost is selected from
#' @param titre_ceiling_gradient (1-A)/B; Where A is the proportion of the full boost received at or above the titre_ceiling_threshold (B)
#' @param titre_ceiling_threshold The maximum titre level an individual can have before their secondary boost is limited in size
#' @param distribution Specifies the distribution with which to randomly sample the boosting parameter
#'
#' @return A boost value gets printed
#' @export
#'
#' @examples simulate_boost(0, 5000, 400,  0.5/200, 200, "log-normal")
simulate_boost <- function(current_titre, boost_mean, boost_sd, titre_ceiling_gradient, titre_ceiling_threshold, distribution){
  boost <- draw_parameters(boost_mean, boost_sd, distribution)
  titre_threshold <- min(current_titre, titre_ceiling_threshold)
  boost <- boost*(1-titre_ceiling_gradient*titre_threshold)
  boost
}
