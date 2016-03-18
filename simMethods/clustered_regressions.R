
################################################################################
# Clustered regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] for r = 1:p[1] is a noisy realization of a much smaller number of
# b[1], ..., b[k]
################################################################################

## ---- libraries ----
library("simData")
library("Matrix")

## ---- opts ----
opts <- list()
opts$n <- 150
opts$p <- c(60, 20) # 60 response dimensions, 20 features
opts$k <- 4
opts$sigma_s0 <- 2
opts$sigma_psi <- .1
opts$sigma_x0 <- 2
opts$sigma_y <- 4

## ---- fixed-params ----
S <- matnorm(opts$p[2], opts$k, opts$sigma_s0)
Psi <- matnorm(opts$p[2], opts$p[2], opts$sigma_psi)
Psi <- t(Psi) %*% Psi # random wishart

## ---- random-params ----
W <- sample(1:opts$k, opts$p[1], replace = T) # each equally likely
dfw <- data.frame(w = as.factor(W))
W <- t(model.matrix(~ w - 1, dfw, intercept = F))

E <- matnorm(opts$p[2], opts$p[1], 1)
B <- S %*% W + sqrtm(Psi) %*% E

W[, 1:8]
pairs(B[, 1:8], B[, 1:8])

## ---- generate-data ----
X <- matnorm(opts$n, opts$p[2], opts$sigma_x0)
Y <- X %*% B + matnorm(opts$n, opts$p[1], opts$sigma_y)
plot(S[, which(W[, 1] == 1)], coef(lm(Y[, 1] ~ X - 1)))
