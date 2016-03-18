
################################################################################
# Generate data from low-rank regression model
################################################################################

## ---- libraries ----
library("simData")
library("corpcor")
library("expm")

## ---- opts ----
opts <- list()
opts$n <- 250
opts$p <- c(60, 20)
opts$k <- 4
opts$sigma_b0 <- 1
opts$sigma_b_noise <- .5
opts$sigma_x0 <- 1
opts$sigma_y <- 4

## ---- make-B ----
U <- matnorm(opts$p[2], opts$k, opts$sigma_b0)
V <- matnorm(opts$p[1], opts$k, opts$sigma_b0)
B <- U %*% t(V)

## ---- make-regression ----
X <- matnorm(opts$n, opts$p[2], opts$sigma_x0)
E <- matnorm(opts$n, opts$p[1], opts$sigma_y)
Y <- X %*% B + E

## ---- make-test-data ----
train_ix <- sample(1:opts$n, 3/4 * opts$n)
Y_test <- Y[-train_ix, ]
X_test <- X[-train_ix, ]
Y <- Y[train_ix, ]
X <- X[train_ix, ]

## ---- naive-regression ----
indep_regs <- lapply(seq_len(ncol(Y)), function(i) lm(Y[, i] ~ X))

plot(X[, 1], Y[, 1])
plot(X[, 2], Y[, 1])
summary(indep_regs[[1]])

## ---- red-rank ----
cancor_res <- cancor(X, Y)
SigmaXX <- cov(X)
SigmaXY <- cov(X, Y)
SigmaYY <- cov(Y)

K <- sqrtm(solve(SigmaXX)) %*% SigmaXY %*% sqrtm(solve(SigmaYY))
V <- svd(K)$v
Pv <- V[, 1:opts$k] %*% t(V[, 1:opts$k])

Q <- qr.Q(qr(X))
R <- qr.R(qr(X))
B_ols <- solve(R) %*% t(Q) %*% Y

Y_hat_ols <- X %*% B_ols
Y_hat_rr <- X %*% B_ols %*% Pv

plot(Y_hat_ols[, 1:2], asp = 1)
points(Y_hat_rr[, 1:2], col = "red")

plot(Y[, 1], Y_hat_ols[, 1], asp = 1, main = "training...")
points(Y[, 1], Y_hat_rr[, 1], col = "red")

## ---- test-red-rank ----
Y_ols_test <- X_test %*% B_ols
Y_rr_test <- X_test %*% B_ols %*% Pv

plot(Y_test[, 1], Y_ols_test[, 1], asp = 1, main = "test")
points(Y_test[, 1], Y_rr_test[, 1], col = "red")

sqrt(mean((Y_ols_test - Y_test)^2))
sqrt(mean((Y_rr_test - Y_test)^2)) # reduced rank does worse?

# an observation: for n in the 150 - 500 range, reduced rank does noticeable
# better than independent regressions. For n larger than this, the two methods
# are practically indistinguishable (there's no need to pool). But
# interestingly, when n is too small, the reduced rank method breaks down,
# because the covariance estimates used to find the canonical correlation axes
# is pretty unreliable.
