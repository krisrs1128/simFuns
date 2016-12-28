
#' @title Return RBF kernel functions
kernel_fun_factory <- function(sigma) {
  function(x, y) {
    (1 / (sqrt(2 * pi * sigma))) * exp(-1 / (2 * sigma) * (x - y) ^ 2)
  }
}

#' @title Generate a Kernel matrix for  GP
#' @export
kernel_matrix <- function(times,  sigma = 1) {
  K <- matrix(0, length(times), length(times))
  f <- kernel_fun_factory(sigma)
  for (i in seq_along(times)) {
    for (j in seq_along(times)) {
      K[i, j] <- f(times[i], times[j])
    }
  }

  K
}

#' @title Take a matrix square root
#' @export
sqrt_mat <- function(X) {
  X_svd <- svd(X)
  X_svd$u %*% diag(sqrt(X_svd$d)) %*% X_svd$v
}

#' @title Generate GP data
#'
#' @examples
#' times <- seq(0, 1, length.out = 100)
#' plot(times, gp_data(times, 0.01, 0.1))
#' @export
gp_data <- function(times, bandwidth = 1, noise = 0) {
  K <- kernel_matrix(times, bandwidth)
  sqrt_mat(K) %*% matnorm(length(times), 1, 1) +
    matnorm(length(times), 1, noise)
}
