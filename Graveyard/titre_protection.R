#' Calculate The Risk Of Infection At A Given Titre
#'
#' @param titre An individual's current titre level
#' @param alpha1 The titre value at which the probability of infection is 0.5
#' @param beta1 This value dictates the shape of the curve; smaller values like 0.001 will produce smoother curves
#'
#' @return The risk of infection is returned
#' @export
#'
#' @examples
titre_protection <- function(titre, alpha1, beta1){
  risk <- 1 - 1/(1 + exp(beta1*(titre - alpha1)))
  risk
}
