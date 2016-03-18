
################################################################################
# Contaminated regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] = eta + eps for r = 1:p[1]. That is, the regression coefficients are
# noisy realizations of some true regression coefficient
################################################################################

## ---- libraries ----
library("simData")
library("Matrix")

## ---- opts ----
opts <- list()
opts$n <- 150
opts$p <- c(60, 20) # 60 response dimensions, 20 features
opts$sigma_eps <- 1
opts$sigma_mu0 <- 3
opts$sigma_x0 <- 1
opts$sigma_psi <- .1
opts$sigma_y <- 1

## ---- fixed-params ----
Psi <- matnorm(opts$p[2], opts$p[2], opts$sigma_psi)
Psi <- Psi %*% t(Psi) # random wishart
mu <- rnorm(opts$p[2], opts$sigma_mu0)

## ---- random-params ----
E <- matnorm(opts$p[2], opts$p[1], 1)
B <- mu %*% t(rep(1, opts$p[1])) + sqrtm(Psi) %*% E
pairs(B[, 1:8], B[, 1:8])

## ---- generate-data ----
X <- matnorm(opts$n, opts$p[2], opts$sigma_x0)
Y <- X %*% B + matnorm(opts$n, opts$p[1], opts$sigma_y)
plot(X[, 1], Y[, 1])
abline(a = B[1, 1], b = 0, col = "red")

coef(lm(Y[, 1] ~ X))
B[, 1]
