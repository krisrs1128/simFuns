
################################################################################
# Comparing methods on data generated from probabilistic canonical correlation
# analysis model
################################################################################

## ---- libraries ----
library("simData")

## ---- opts ----
opts <- list()
opts$n <- 150
opts$p <- c(20, 60)
opts$k_shared <- 4
opts$k_X <- 1
opts$k_Y <- 1
opts$sigma_w0 <- 1
opts$sigma_s0 <- 2 # same across all tables
opts$sigma_x <- .5
opts$sigma_y <- .5

## ---- fixed-params ----
S <- list()
S[["X_shared"]] <- matnorm(opts$p[2], opts$k_shared, opts$sigma_s0)
S[["Y_shared"]] <- matnorm(opts$p[1], opts$k_shared, opts$sigma_s0)
S[["X_unique"]] <- matnorm(opts$p[2], opts$k_X, opts$sigma_s0)
S[["Y_unique"]] <- matnorm(opts$p[1], opts$k_Y, opts$sigma_s0)

## ---- random-params ----
W <- list()
W[["shared"]] <- matnorm(opts$n, opts$k_shared, opts$sigma_w0)
W[["X"]] <- matnorm(opts$n, opts$k_X, opts$sigma_w0)
W[["Y"]] <- matnorm(opts$n, opts$k_Y, opts$sigma_w0)

## ---- generate-data ----
E <- list()
E[["X"]] <- matnorm(opts$n, opts$p[2], opts$sigma_x)
E[["Y"]] <- matnorm(opts$n, opts$p[2], opts$sigma_y)

Y <- W[["shared"]] %*% t(S[["Y_shared"]]) + W[["Y"]] %*% t(S[["Y_unique"]])
X <- W[["shared"]] %*% t(S[["X_shared"]]) + W[["X"]] %*% t(S[["X_unique"]])

pairs(cbind(Y[, 1:3], X[, 1:3]))
