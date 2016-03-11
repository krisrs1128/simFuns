
################################################################################
# Simulate data according to Daniela Witten's PMD paper
################################################################################

## ---- packages ----
library("reshape2")
library("ggplot2")
library("dplyr")
library("PMA")
theme_set(theme_bw())

## ---- opts ----
opts <- list()
opts$n <- 504 # want divisible by 8
opts$p <- 20
opts$k <- 2
opts$l <- 2
opts$sigma0 <- 5
opts$sigma <- 1

## ---- source-scores ----
S <- qr.Q(qr(matnorm(opts$p, opts$k, opts$sigma0))) # same source between the two matrices
W <- replicate(opts$l, matrix(0, opts$n, opts$k), simplify = F)
W[[1]][, 1] <- c(rep(10, opts$n / 8), rep(-10, opts$n / 8),
                  rep(10, opts$n / 8), rep(-10, opts$n / 8),
                  rep(0, opts$n / 2))
W[[1]][, 2] <- c(rep(10, opts$n / 4), rep(-10, opts$n / 4),
                  rep(0, opts$n / 2))
W[[2]][, 1] <- c(rep(0, opts$n / 4), rep(10, opts$n / 2), rep(-10, opts$n / 4))
W[[2]][, 2] <- c(rep(-10, opts$n / 4), rep(10, opts$n / 4), rep(-10, opts$n / 4), rep(10, opts$n / 4))

## ---- generate-data ----
X <- common_source_model(W, S, opts)
X <- lapply(X, scale)
pairs(X[[1]][, 1:4])

## ---- separate-pcas ----
pca_sep <- lapply(X, princomp)
pca_sep_scores <- lapply(pca_sep, function(x) x$scores)

D <- melt(list(Mu = Mu, pca_sep = pca_sep_scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- concatenated-pca ----
pca_concat <- princomp(do.call(cbind, X))

D <- melt(list(Mu = Mu, pca_concat = pca_concat$scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")
ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_wrap(type ~ table ~ comp, scale = "free_y")

## ---- pmd-unordered ----
pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 10, ncomponents = 3)
D2 <- melt(list(Mu = Mu, pmd_sep = pmd_res$ws))
colnames(D2) <- c("ix", "comp", "value", "table", "type")
ggplot(D2 %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- pmd-ordered ----
pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 1, type = "ordered",
                    ncomponents = 3)
D3 <- melt(list(Mu = Mu, pmd_sep = pmd_res$ws))
colnames(D3) <- c("ix", "comp", "value", "table", "type")
ggplot(D3 %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- mfa ----
library("FactoMineR")
mfa_res <- MFA(do.call(cbind, X), group = sapply(X, ncol))
D <- melt(list(Mu = Mu,
               mfa_group_one = mfa_res$separate.analyses$group.1$ind$coord[, 1:3],
               mfa_group_two = mfa_res$separate.analyses$group.2$ind$coord[, 1:3]))
str(mfa_res)
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")
