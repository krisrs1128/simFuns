
pdf_factory <- function(A) {
  function(x) {
    A %*% x + rnorm(length(x))
  }
}

#' @title Generic simulation from a dynamical system
#' @example
#' u <- dynamical_system(pdf_factory(diag(2)), c(0, 0))
dynamical_system <- function(pdf, x_init = 0, len = 100) {
  x <- matrix(nrow = length(x_init), ncol = len)
  for (i in seq_len(len)) {
    if (i == 1) {
      x[, i] <- pdf(x_init)
    } else {
      x[, i] <- pdf(x[, i - 1])
    }
  }
  x
}

ssm_sample <- function(pdf_u, u_init = 0, B = NULL, R = NULL, len = 100) {
  u <- dynamical_system(pdf_u, u_init, len)
  K <- nrow(u)

  if (is.null(B)) {
    B <- diag(K)
  }
  p <- nrow(B)

  if (is.null(R)) {
    R <- diag(p)
  }

  B %*% u + sqrt_mat(R) %*% matnorm(p, len, 1)
}

#' @title Generate many samples from a state space model
#' @description This doesn't allow any input data (yet)
#' @examples
#' pdf_u <- pdf_factory(matnorm(2, 2, 1))
#' B <- matnorm(4, 2, .1)
#' ssm_data(10, pdf_u, B, u_init = rep(0, 2))
ssm_data <- function(n, pdf, u_init = 0, B = NULL, R = NULL, len = 100) {
  replicate(n, ssm_sample(pdf, u_init, B, R, len))
}
