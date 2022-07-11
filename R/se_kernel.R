#' Kernel
#'
#' @param x
#' @param y
#' @param sigma
#' @param length
#'
#' @return
#' @export
#'
#' @examples
se_kernel <- function(x, y, sigma = 1, length = 1) {
  sigma^2 * exp(- (x - y)^2 / (2 * length^2))
}
