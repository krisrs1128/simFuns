
################################################################################
#
################################################################################

merge_sparse_ppls_opts <- function(opts = list()) {
  opts <- merge_ppls_opts(opts)
  if(is.null(opts$zero_mask)) {
    opts$zero_mask <- list(
      "u" = replicate(opts$l, matrix(1, opts$n, opts$k_shared), simplify = F),
      "w" = replicate(opts$l, matrix(1, opts$n, opts$k_unique), simplify = F),
      "v" = replicate(opts$l, matrix(1, opts$p, opts$k_shared), simplify = F),
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

  V_sparse <- sparsify(ppls_res$params$V, opts$zero_mask_v)
  for(l in seq_along(ppls_res$params$U)) {
    ppls_res$params$U[[l]] <- sparsify(ppls_res$params$U[[l]], opts$zero_masks[["u"]])
    ppls_res$params$V[[l]] <- sparsify(ppls_res$params$V[[l]], opts$zero_masks[["v"]])
    ppls_res$params$S[[l]] <- sparsify(ppls_res$params$S[[l]], opts$zero_masks[["S"]])
    X[[l]] <- params$W[[l]] %*% t(params$S[[l]]) +
      params$U[[l]] %*% t(params$V) +
      matnorm(opts$n, opts$p, opts$sigma)
  }
  list(X = X, params = ppls_res$params)
}
