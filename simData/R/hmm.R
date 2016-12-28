
rmultinoulli <- function(p, n = 1) {
  u <- runif(n)
  sum_p <- cumsum(p)
  stopifnot(sum_p[length(p)] > .9999 && sum_p[length(p)] < 1.0001)

  z <- rep(1, n)
  for (i in length(p):1) {
    z[u <= sum_p[i]] <- i
  }
  z
}

markov_chain <- function(P, pi_init = NULL, length = 100) {
  z <- vector(length = length)
  K <- nrow(P)

  if (is.null(pi_init)) {
    pi_init <- rep(1 / K, K)
  }

  z[1] <- rmultinoulli(pi_init)

  for (i in 2:length) {
    z[i] <- rmultinoulli(P[z[i - 1], ])
  }

  z
}

#' @title Generate many markov chains
#'
#' @examples
#' # pretty diagonally dominant
#' P <- matrix(runif(4 * 4), 4, 4) + 10 * diag(4)
#' P <- diag(1 / rowSums(P)) %*% P
#' matplot(t(markov_chains(P, n = 3)))
#' @export
markov_chains <- function(P, pi_init = NULL, n = 1, length = 100) {
  t(replicate(n,  markov_chain(P, pi_init, length)))
}

#' @title Generate data from an HMM
#' @examples
#' P <- matrix(runif(4 * 4), 4, 4) + 10 * diag(4)
#' P <- diag(1 / rowSums(P)) %*% P
#' obs_densities <- lapply(1:4, function(sigma) { function(n) { rnorm(n, 0, sigma) }})
#' names(obs_densities) <- 1:4
#' x <- hmm_data(P, obs_densities, 4)
#' matplot(t(x))
#' @export
hmm_data <- function(P, obs_densities, n = 1, pi_init = NULL, length = 100) {
  z <- markov_chains(P, pi_init, n, length)
  x <- matrix(0, nrow(z), ncol(z))
  K <- nrow(P)

  for (k in seq_len(K)) {
    n_obs <- sum(z == k)
    x[z == k] <- obs_densities[[as.character(k)]](n_obs)
  }
  x
}
