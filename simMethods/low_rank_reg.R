
################################################################################
# Generate data from low-rank regression model
################################################################################

## ---- libraries ----
library("simData")

## ---- opts ----
opts <- list()
opts$n <- 100
opts$p <- c(40, 20)
opts$k <- 4
opts$sigma_b0 <- 1
opts$sigma_b_noise <- .2
opts$sigma_x0 <- 1
opts$sigma_y <- 1

## ---- make-B ----
U <- matnorm(opts$p[2], opts$k, opts$sigma_b0)
V <- matnorm(opts$p[1], opts$k, opts$sigma_b0)
B <- U %*% t(V)

## ---- make-regression ----
X <- matnorm(opts$n, opts$p[2], opts$sigma_x0)
E <- matnorm(opts$n, opts$p[1], opts$sigma_y)
Y <- X %*% B + E

## ---- naive-regression ----
indep_regs <- lapply(seq_len(ncol(Y)), function(i) lm(Y[, i] ~ X))

plot(X[, 1], Y[, 1])
plot(X[, 2], Y[, 1])
summary(indep_regs[[1]])
