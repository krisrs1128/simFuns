
################################################################################
# Simulate data according to Daniela Witten's PMD paper
################################################################################

## ---- sparse-sim-packages ----
# List of packages for session
.packages = c("reshape2",
              "ggplot2",
              "plyr",
              "dplyr", 
              "simData",
              "PMA", 
              "vegan",
              "FactoMineR")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.packages[!.inst], repos='http://cran.rstudio.com/')
}

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
theme_set(small_theme())
set.seed(04032016)

cat("\014")  # Clear console

rm(list=ls()) # Delete all existing variables
graphics.off() # Close all open plots

## ---- sparse-sim-opts ----
opts <- list()
opts$n <- 504 # want divisible by 8
opts$p <- 20
opts$k <- 2
opts$l <- 2
opts$sigma0 <- 5
opts$sigma <- 1

## ---- sparse-sim-source-weights ----
S <- qr.Q(qr(matnorm(opts$p, opts$k, opts$sigma0))) # same source between the two matrices
W <- replicate(opts$l, matrix(0, opts$n, opts$k), simplify = F)
W[[1]][, 1] <- c(rep(10, opts$n / 8), rep(-10, opts$n / 8),
                  rep(10, opts$n / 8), rep(-10, opts$n / 8),
                  rep(0, opts$n / 2))
W[[1]][, 2] <- c(rep(10, opts$n / 4), rep(-10, opts$n / 4),
                  rep(0, opts$n / 2))
W[[2]][, 1] <- c(rep(0, opts$n / 4), rep(10, opts$n / 2), rep(-10, opts$n / 4))
W[[2]][, 2] <- c(rep(-10, opts$n / 4), rep(10, opts$n / 4), rep(-10, opts$n / 4), rep(10, opts$n / 4))

## ---- sparse-sim-plot-weights-setup ----
mW <- melt(W)
colnames(mW) <- c("i", "k", "w", "table")
mW$k <- paste0("Latent Dimension ", mW$k)
mW$table[mW$table == 1] <- "X"
mW$table[mW$table == 2] <- "Y"

## ---- sparse-sim-plot-weights ----
ggplot(mW) +
  geom_point(aes(x = i, y = w)) +
  xlab("sample index") + 
  facet_grid(table ~ k) +
  ggtitle(paste0("True latent weights W"))

## ---- sparse-sim-generate-data ----
X <- common_source_model(W, S, opts)

## ---- sparse-sim-plot-data-x ----
colnames(X[[1]]) <- paste0("X", seq_len(ncol(X[[1]])))
colnames(X[[2]]) <- paste0("Y", seq_len(ncol(X[[2]])))
x_ix <- sample(seq_len(ncol(X[[1]])), 4)
y_ix <- sample(seq_len(ncol(X[[1]])), 4)
pairs(X[[1]][, x_ix], asp = 1, main = "Four columns of X")

## ---- sparse-sim-plot-data-y ----
pairs(X[[2]][, y_ix], asp = 1, main = "Four columns of Y")

## ---- sparse-sim-plot-data-xy ----
pairs(cbind(X[[1]][, x_ix[1:2]], X[[2]][, y_ix[1:2]]), asp = 1,
      main = "Two columns of X vs. Two columns of Y")

## ---- sparse-sim-pmd-unordered ----
X <- lapply(X, scale)
pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 10, ncomponents = 3)

## ---- sparse-sim-pmd-unordered-plot ----
mW_hat <- melt(pmd_res$ws)
colnames(mW_hat) <- c("i", "k", "w", "table")
mW_hat <- mW_hat %>% filter(k < 4)
mW_hat$k <- paste0("Recovered Dimension ", mW_hat$k)
mW_hat$table[mW_hat$table == 1] <- "X"
mW_hat$table[mW_hat$table == 2] <- "Y"

ggplot(mW_hat) +
  geom_point(aes(x = i, y = w), alpha = 0.6, size = 1) +
  facet_grid(table ~ k) +
  ggtitle(expression(paste("Recovered Weights ", hat(W))))

## ---- sparse-sim-pmd-ordered ----
pmd_res <- MultiCCA(lapply(X, function(x) t(x)), penalty = 1, type = "ordered",
                    ncomponents = 3)
D3 <- melt(list(Mu = Mu, pmd_sep = pmd_res$ws))
colnames(D3) <- c("ix", "comp", "value", "table", "type")
ggplot(D3 %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- sparse-sim-separate-pcas ----
pca_sep <- lapply(X, princomp)
pca_sep_scores <- lapply(pca_sep, function(x) x$scores)
Mu <- W

D <- melt(list(Mu = Mu, pca_sep = pca_sep_scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- sparse-sim-concatenated-pca ----
pca_concat <- princomp(do.call(cbind, X))

D <- melt(list(Mu = Mu, pca_concat = pca_concat$scores))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6, size = 1) +
  facet_wrap(type ~ comp, scale = "free_y")

## ---- sparse-sim-mfa ----
mfa_res <- MFA(do.call(cbind, X), group = sapply(X, ncol))
D <- melt(list(Mu = Mu,
               mfa_group_one = mfa_res$separate.analyses$group.1$ind$coord[, 1:3],
               mfa_group_two = mfa_res$separate.analyses$group.2$ind$coord[, 1:3]))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")

## ---- sparse-sim-cca ----
cca_res <- vegan::CCorA(X[[1]], X[[2]])
D <- melt(list(Mu = Mu,
               cca_x1 = cca_res$Cy[, 1:3],
               cca_x2 = cca_res$Cx[, 1:3]))
colnames(D) <- c("ix", "comp", "value", "table", "type")

ggplot(D %>% filter(comp < 4)) +
  geom_point(aes(x = ix, y = value, col = as.factor(comp), shape = type), alpha = 0.6,
             size = 1) +
  facet_grid(type ~ table ~ comp, scale = "free_y")
