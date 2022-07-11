#' Given X Coordinates, Take N Draws From Kernel Function At Those Points
#'
#' @param x The number of coordinates
#' @param N  The number of draws to take
#' @param seed
#' @param kernel_fn
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
draw_samples <- function(x, N, seed = 1, kernel_fn, ...) {
  Y <- matrix(NA, nrow = length(x), ncol = N)
  #set.seed(seed)
  for (n in 1:N) {
    K <- cov_matrix(x, kernel_fn, ...)
    Y[, n] <- MASS::mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
  }
  Y
}
