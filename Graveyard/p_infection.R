#' Generate The Probability Of Infection At A Given Titre
#'
#' @param phi The probability of infection
#' @param titre An individual's current titre level
#' @param alpha1 The titre value at which the probability of infection is 0.5
#' @param beta1 This value dictates the shape of the curve; smaller values like 0.001 will produce smoother curves
#'
#' @return The probability of infection is returned
#' @export
#'
#' @examples
p_infection <- function(phi, titre, alpha1, beta1){
  p <- phi*(1-titre_protection(titre, alpha1 , beta1))
  p
}
