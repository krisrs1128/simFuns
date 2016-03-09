
################################################################################
# Generate data according to probabilistic PLS model
################################################################################

#' @title Default options for PPLS simulation
#' @param opts  A (potentially partially specified) list of options to use when
#' simulating from the probabilistic PLS model. The potential arguments are \cr
#'    $k_shared The latent dimension for the shared sources. \cr
#'    $k_unique The latent dimensions for the table-unique sources.
#'    $n The total number of samples, for each table.
#'    $p The number of features, for each table.
#'    $sigma The noise variance.
#'    $sigma0 The variance when generating scores / factors.
#' @return opts The original opts with defaults filled in.
#' @export
merge_ppls_opts <- function(opts) {
  default_opts <- list()
  default_opts$k_shared <- 2
  default_opts$k_unique <- c(2, 2)
  default_opts$n <- 100
  default_opts$p <- 10
  default_opts$sigma <- 1
  default_opts$sigma0 <- 4
  modifyList(default_opts, opts)
}

#' @title Helper to create random normal matrices
#' @export
matnorm <- function(n, p, sigma) {
  matrix(rnorm(n * p, 0, sd = sigma), n, p)
}

#' @title Probabilistic PLS simulation
#' @description Generate L tables X[[l]], so that
#'      X[[l]] = W[[l]] S[[l]]' + U * V' + E[[l]]
#' where the W's and S's are k[[l]]-dimensional latent scores / factors,
#' U and V are shared k_shared-dimensional scores / factors, generated from
#' spherical normals with variance opts$sigma times the identity, and E is
#' spherical normal with variance opts$sigma times the identity.
#'
#' Note this is not the "true" ppls model -- there the prior is only on the
#' latent factors / sources, not scores / coordinates. Here we just need a way
#' to generate those coordinates, so we use a normal.
#'
#' @param opts Options for the PPLS simulation. See merge_ppls_opts().
#' @examples
#' ppls_res <- ppls_shared()
#' X <- do.call(cbind, ppls_res$X)
#' res <- princomp(X)
#' plot(res$scores[, 1], ppls_res$params$U[, 1])
#' @export
ppls_shared <- function(opts = list()) {
  opts <- merge_ppls_opts(opts)
  params <- list(W = list(), S = list(), U = NULL, V = NULL)
  X <- list()

  # data unique to each table
  params$U <- matnorm(opts$n, opts$k_shared, opts$sigma0)
  params$V <- matnorm(opts$p, opts$k_shared, opts$sigma0)
  for(l in seq_along(opts$k_unique)) {
    params$W[[l]] <- matnorm(opts$n, opts$k_unique[l], opts$sigma0)
    params$S[[l]] <- matnorm(opts$p, opts$k_unique[l], opts$sigma0)
    X[[l]] <- params$W[[l]] %*% t(params$S[[l]]) +
      params$U %*% t(params$V) +
      matnorm(opts$n, opts$p, opts$sigma)
  }

  list(X = X, params = params)
}
