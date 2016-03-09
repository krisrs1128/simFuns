
################################################################################
# Simulations for classical multivariate methods
################################################################################

## ---- packages ----
library("simData")
library("ggplot2")
library("dplyr")
library("reshape2")
theme_set(theme_bw())

## ---- load-data ----
data(pls_two_table)
data(pls_four_table)

## ---- pca ----
X_merged2 <- do.call(cbind, pls_two_table$X)
X_merged4 <- do.call(cbind, pls_four_table$X)

pls_svd2 <- svd(X_merged2)
pls_svd4 <- svd(X_merged4)

UD_pls_2 <- pls_svd2$u %*% diag(pls_svd2$d)
UD_pls_4 <- pls_svd4$u %*% diag(pls_svd4$d)

V_pls_2 <- pls_svd2$v
V_pls_4 <- pls_svd4$v

## ---- plot-true-scores ----
true_scores_2 <- melt(pls_two_table$params$W, varnames = c("index", "dimension")) %>%
  dcast(index + L1 ~ dimension)
colnames(true_scores_2) <- c("index", "table", "dim1", "dim2")

ggplot(true_scores_2) +
  geom_point(aes_string(x = "dim1", y = "dim2")) +
  facet_wrap(~table) +
  ggtitle("True PLS scores")

true_scores_4 <- melt(pls_four_table$params$W, varnames = c("index", "dimension")) %>%
  dcast(index + L1 ~ dimension)
colnames(true_scores_4) <- c("index", "table", "dim1", "dim2")

ggplot(true_scores_4) +
  geom_point(aes_string(x = "dim1", y = "dim2")) +
  facet_wrap(~table) +
  ggtitle("True PLS scores")

## ---- compare-scores ----
dim(UD_pls_2)
head(cor(UD_pls_2, pls_two_table$params$W[[1]]))
plot(UD_pls_2[, 4], pls_two_table$params$W[[1]][, 1])
