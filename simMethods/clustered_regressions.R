
################################################################################
# Clustered regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] for r = 1:p[1] is a noisy realization of a much smaller number of
# b[1], ..., b[k]
################################################################################

## ---- libraries ----
library("simData")
library("expm")
library("plyr")
library("dplyr")

## ---- opts ----
opts <- list()
opts$n <- 1000
opts$p <- c(100, 2) # 60 response dimensions, 20 features
opts$k <- 3
opts$sigma_x0 <- 2
opts$sigma_b <- .8
opts$sigma_y <- 2

## ---- fixed-params ----
S <- matrix(c(0, 2, 3, 1, 2, 4), nrow = 2)
Psi <- diag(opts$p[2])

## ---- random-params ----
W <- sample(1:opts$k, opts$p[1], replace = T) # each equally likely
dfw <- data.frame(w = as.factor(W))
W <- t(model.matrix(~ w - 1, dfw, intercept = F))

E <- matnorm(opts$p[2], opts$p[1], opts$sigma_b)
B <- S %*% W + sqrtm(Psi) %*% E

plot(t(B), asp = 1)

## ---- generate-data ----
x_list <- list()
y_list <- list()
for(r in seq_len(opts$p[1])) {
  x_list[[r]] <- matnorm(opts$n, opts$p[2], opts$sigma_x0)
  y_list[[r]] <- x_list[[r]] %*% B[, r] + rnorm(opts$n, 0, opts$sigma_y)
}

plot(x_list[[1]][, 1], y_list[[1]])

phi <- rep(1 / opts$k, opts$k)
sigma <- opts$sigma_y # assume we know sigma for now

S0 <- matrix(c(1.5, 2.5, 2.5, 1.5, .5, 2.0), nrow = 2)
Psi0 <- diag(runif(opts$p[2], .8, 1.2))
data_list <- list(x_list = x_list, y_list = y_list)
param_list <- list(S = S0, Psi = Psi0, phi = phi, sigma = sigma)

vb_res <- vb_multinom(data_list, param_list, 3)

pi_guess <- t(vb_res$var_list$pi)
plot(t(W)[, 1], pi_guess[, 1])
plot(t(W)[, 2], pi_guess[, 2])
m <- t(vb_res$var_list$m)

plot(m)
clust <- apply(pi_guess, 1, which.max)
points(m[clust == 1, ], col = 'purple')
points(m[clust == 2, ],  col = 'blue')

plot(B[1, ], vb_res$var_list$m[1, ], asp = 1)
plot(B[2, ], vb_res$var_list$m[2, ], asp = 1)

plot(vb_res$elbo)

plot(295:320, vb_res$elbo[295:320])

plot(log(-vb_res$elbo))
#abline(h = 0)
min(diff(vb_res$elbo))
vb_res$param_list$Psi

table(apply(W, 2, which.max), clust)

vb_res$var_list$m[, 1]
plot(x_list[[1]][, 1], y_list[[1]])
abline(a = 0, b = B[1, 1])
B1_hat <- vb_res$param_list$S %*% vb_res$var_list$pi[, 1]
abline(a = 0, b = B1_hat[1], col = "red")

