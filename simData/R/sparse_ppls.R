
################################################################################
#
################################################################################

merge_sparse_ppls_opts <- function(opts = list()) {
  opts <- merge_ppls_opts(opts)
  if(is.null(opts$zero_mask)) {
    opts$zero_mask <- list(
      "v" = matrix(1, opts$p, opts$k_shared),
      "u" = replicate(opts$l, matrix(1, opts$n, opts$k_shared), simplify = F),
      "w" = replicate(opts$l, matrix(1, opts$n, opts$k_unique), simplify = F),
      "s" = replicate(opts$l, matrix(1, opts$p, opts$k_unique), simplify = F)
    )
  }
  opts
}

sparsify <- function(X, zero_mask) {
  stopifnot(dim(X) == dim(zero_mask))
  X[zero_mask == 0] <- 0
  X
}

# zero_mask_u = list of matrices same dimension as u, with 0's and 1's
# indicating whether to zero out certain entries
sparse_ppls_shared <- function(opts = list()) {
  opts <- merge_sparse_ppls_opts(opts)
  ppls_res <- ppls_shared(opts)

  ppls_res$params$V <- sparsify(ppls_res$params$V, opts$zero_mask[["v"]])
  for(l in seq_along(ppls_res$params$U)) {
    ppls_res$params$U[[l]] <- sparsify(ppls_res$params$U[[l]], opts$zero_masks[["u"]][[l]])
    ppls_res$params$W[[l]] <- sparsify(ppls_res$params$W[[l]], opts$zero_masks[["W"]][[l]])
    ppls_res$params$S[[l]] <- sparsify(ppls_res$params$S[[l]], opts$zero_masks[["S"]][[l]])
    X[[l]] <- ppls_res$params$W[[l]] %*% t(ppls_res$params$S[[l]]) +
      ppls_res$params$U[[l]] %*% t(ppls_res$params$V) +
      matnorm(opts$n, opts$p, opts$sigma)
  }
  list(X = X, params = ppls_res$params)
}
