
################################################################################
# Simulations for classical multivariate methods
################################################################################

## ---- packages ----
library("simData")

## ---- load-data ----
data(pls_two_table)
data(pls_four_table)

## ---- pca ----
X_merged2 <- do.call(cbind, pls_two_table$X)
X_merged4 <- do.call(cbind, pls_four_table$X)

pls_svd2 <- svd(X_merged2)
pls_svd4 <- svd(X_merged4)

UD_pls_2 <- pls_svd2$u[, 1:2] %*% diag(pls_svd2$d[1:2])
UD_pls_4 <- pls_svd4$u[, 1:2] %*% diag(pls_svd4$d[1:2])

V_pls_2 <- pls_svd2$v[, 1:2]
V_pls_4 <- pls_svd4$v[, 1:2]
