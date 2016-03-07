
################################################################################
# Generate data according to probabilistic PCA model
################################################################################

#' @title Default options for PPCA simulation
#' @param opts  A (potentially partially specified) list of options to use when
#' simulating from the probabilistic PCA model. The potential arguments are \cr
#'    $k_shared The latent dimension for the shared sources. \cr
#'    $k_unique The latent dimensions for the table-unique sources.
#'    $n The total number of samples, for each table.
#'    $p The number of features, for each table.
#'    $sigma The noise variance.
#'    $sigma0 The variance when generating scores / factors.
#' @return opts The original opts with defaults filled in.
#' @export
merge_ppca_opts <- function(opts) {
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

#' @title Probabilistic PCA simulation
#' @description Generate L tables X[[l]], so that
#'      X[[l]] = W[[l]] S[[l]]' + U * V' + E[[l]]
#' where the W's and S's are k[[l]]-dimensional latent scores / factors,
#' U and V are shared k_shared-dimensional scores / factors, generated from
#' spherical normals with variance opts$sigma times the identity, and E is
#' spherical normal with variance opts$sigma times the identity.
#'
#' Note this is not the "true" ppca model -- there the prior is only on the
#' latent factors / sources, not scores / coordinates. Here we just need a way
#' to generate those coordinates, so we use a normal.
#'
#' @param opts Options for the PPCA simulation. See merge_ppca_opts().
#' @examples
#' ppca_res <- ppca_shared()
#' X <- do.call(cbind, ppca_res$X)
#' res <- princomp(X)
#' plot(res$scores[, 1], ppca_res$params$U[, 1])
#' @export
ppca_shared <- function(opts = list()) {
  opts <- merge_ppca_opts(opts)
  params <- list(W = list(), S = list(), U = NULL, V = NULL)
  X <- list()

  # data unique to each table
  params$U <- matnorm(opts$n, opts$k_shared, opts$sigma0)
  params$V <- matnorm(opts$p, opts$k_shared, opts$sigma0)
  for(l in seq_along(opts$p)) {
    params$W[[l]] <- matnorm(opts$n, opts$k_unique[l], opts$sigma0)
    params$S[[l]] <- matnorm(opts$p, opts$k_unique[l], opts$sigma0)
    X[[l]] <- params$W[[l]] %*% t(params$S[[l]]) +
      matnorm(opts$n, opts$p, opts$sigma)
  }

  list(X = X, params = params)
}

