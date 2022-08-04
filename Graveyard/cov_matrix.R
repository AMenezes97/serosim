#
#' Generate Covariance Matrix For Points In `x` Using Given Kernel Function
#'
#' @param x
#' @param kernel_fn
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cov_matrix <- function(x, kernel_fn, ...) {
  outer(x, x, function(a, b) kernel_fn(a, b, ...))
}
