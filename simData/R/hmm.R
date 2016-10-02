
rmultinoulli <- function(p, n = 1) {
  u <- runif(n)
  sum_p <- cumsum(p)
  stopifnot(sum_p[length(p)] > .999 && sum_p[length(p)] < 1.001)

  z <- rep(1, n)
  for (i in length(p):2) {
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

#' @ generate many markov chains
#'
#' @example
#' # pretty diagonally dominant
#' P <- matrix(runif(4 * 4), 4, 4) + 10 * diag(4)
#' P <- diag(1 / rowSums(P)) %*% P
#' matplot(t(markov_chains(10, P)))
markov_chains <- function(n, P, pi_init = NULL, length = 100) {
  t(replicate(n,  markov_chain(P, pi_init, length)))
}
