
################################################################################
# What happens to ordination methods when we change the difference in number of
# columns across tables?
################################################################################

## ---- libraries ----
library("simData")
library("reshape2")
library("ggplot2")
library("plyr")
library("FactoMineR")
library("dplyr")
theme_set(theme_bw())

## ---- opts ----
opts <- list()
opts$n <- 1000
opts$p <- c(5, 900)
opts$k <- 2
opts$k_unique <- c(2, 5)
opts$l <- 2
opts$sigma0 <- 5
opts$sigma <- 1

## ---- utils ----
table_balance_data <- function(opts) {
  S <- list()
  W <- list()

  for(l in seq_len(opts$l)) {
    S[[l]] <- qr.Q(qr(matnorm(opts$p[[l]], opts$k_unique[[l]], 1)))
    W[[l]] <- matnorm(opts$n, opts$k_unique[[l]], 1)
  }

  W[[1]][, 1] <- rnorm(opts$n, c(rep(-3, opts$n / 2), rep(3, opts$n / 2)))
  W[[2]][, 1] <- rnorm(opts$n, c(rep(3, 3 * opts$n / 4), rep(-3, opts$n / 4)))
  W[[2]][, 2] <- rnorm(opts$n, c(rep(3, opts$n / 4), rep(-3, 3 * opts$n / 4)))

  X <- list()
  for(l in seq_len(opts$l)) {
    X[[l]] <- W[[l]] %*% t(S[[l]]) + matnorm(opts$n, opts$p[[l]], opts$sigma)
  }
  list(W = W, S = S, X = X)
}

## ---- source-scores ----
tb_data <- table_balance_data(opts)

pairs(tb_data$W[[1]])
pairs(tb_data$W[[2]])
pairs(tb_data$X[[1]])

## ---- pca-sep ----
pca_sep <- lapply(tb_data$X, princomp)
pca_sep_scores <- lapply(pca_sep, function(x) x$scores)

D <- melt(list(W = tb_data$W, pca_sep = pca_sep_scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")
ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- pca-concat ----
pca_concat <- princomp(do.call(cbind, tb_data$X))

D <- melt(list(W = tb_data$W, pca_concat = pca_concat$scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")
ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- mfa ----
mfa_res <- MFA(do.call(cbind, tb_data$X), group = sapply(tb_data$X, ncol), graph = F)
D <- melt(list(W = tb_data$W,
               mfa_group_one = mfa_res$separate.analyses$group.1$ind$coord[, 1:3],
               mfa_group_two = mfa_res$separate.analyses$group.2$ind$coord[, 1:3]))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")
