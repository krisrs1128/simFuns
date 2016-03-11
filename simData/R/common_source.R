
################################################################################
# Simple common source model
################################################################################

#' @title Model with common sources, but different scores for each
#' @param W A list of matrices containing scores for each table.
#' @param S A single source matrix underlying all the tables.
#' @return A list of matrices, each with the same sources S but different
#' scores W.
#' @export
common_source_model <- function(W, S, opts) {
  X <- list()
  n <- nrow(W[[1]])
  p <- nrow(S)

  for(l in seq_along(W)) {
    X[[l]] <- W[[l]] %*% t(S) + matnorm(n, p, opts$sigma)
  }

  X
}
