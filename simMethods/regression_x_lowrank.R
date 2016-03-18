
################################################################################
# A simulation Y ~ XB but where X is low rank. This is considered an ideal
# situation for PCA-IV, but PLS may also perform well here.
################################################################################

## ---- libraries ----
library("simData")
library("ggplot2")
library("ade4")
library("pls")
library("reshape2")
library("plyr")
library("dplyr")
theme_set(theme_bw())

## ---- opts ----
opts <- list()
opts$n <- 150
opts$p <- c(15, 30)
opts$k <- 4
opts$sigma <- 1
opts$sigma_x_noise <- sqrt(2 * log(opts$n))
opts$sigma_b <- 1
opts$sigma_y <- 10

## ---- simulate-x ----
U <- matnorm(opts$n, opts$k, opts$sigma)
V <- matnorm(opts$p[2], opts$k, opts$sigma)

X <- U %*% t(V) + matnorm(opts$n, opts$p[2], opts$sigma_x_noise)
plot(svd(X)$d)

## ---- simulate-regression ----
B <- matnorm(opts$p[2], opts$p[1], opts$sigma)
E <- matnorm(opts$n, opts$p[1], opts$sigma_y)
Y <- X %*% B + E
train_ix <- seq_len(3 * opts$n / 4)

X_test <- X[-train_ix, ]
X <- X[train_ix, ]
Y_test <- Y[-train_ix, ]
Y <- Y[train_ix, ]

## ---- naive-regression ----
summary(lm(Y[, 1] ~ X)) # don't want to regress when low rank!

## ---- pca-iv-ade4 ----
run_pca_iv <- function(X, Y, nfX, nfY) {
  pca_Y <- dudi.pca(Y, scan = F, nf = nfY)
  pca_iv_YX <- pcaiv(pca_Y, X, scan = F, nf = nfX)
  print(summary(pca_iv_YX))
  plot(pca_iv_YX)
  pca_iv_YX
}

pcaiv_res <- run_pca_iv(X, Y, opts$k, opts$k)

## ---- pca-iv-by-hand ----
p <- sum(opts$p)
Sigma <- cov(cbind(scale(X), scale(Y)))
SigmaYY <- Sigma[1:opts$p[1], 1:opts$p[1]]
SigmaXX <- Sigma[(opts$p[1] + 1):p, (opts$p[1] + 1):p]
SigmaXY <- Sigma[(opts$p[1] + 1):p, 1:opts$p[1]]

pcaiv_res_direct <- svd(solve(SigmaXX) %*% SigmaXY %*% t(SigmaXY))
V <- pcaiv_res_direct$v
lambda <- pcaiv_res_direct$d

## ---- compare-res ----
plot(lambda)
points(pcaiv_res$eig, col = "red")

## ---- pls ----
pls_res <- plsr(Y ~ X)
plot(pls_res)
coefplot(pls_res)
corrplot(pls_res)

Y_test_hat <- predict(pls_res, X_test)
colnames(Y_test_hat) <- NULL

mY_test <- rbind(data.frame(type = "truth", Y = Y_test, ix = seq_len(nrow(Y_test))),
                 data.frame(type = "fitted", Y = Y_test_hat[, , 30], ix = seq_len(nrow(Y_test)))) %>%
                   melt(id.vars = c("type", "ix")) %>%
                   dcast(ix + variable ~ type)

ggplot(mY_test) +
  geom_text(aes(x = truth, y = fitted, label = ix), size = 3) +
  facet_wrap(~variable)

sqrt(mean((Y_test - Y_test_hat[, , 30]) ^ 2))
sqrt(mean((Y_test - mean(Y_test))^2)) # baseline


