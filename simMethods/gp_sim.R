
################################################################################
# Simulate data according to clustered gaussian processes
################################################################################

## ---- packages ----
library("ggplot2")
library("simData")
library("plyr")
library("dplyr")
library("reshape2")
library("stringr")
theme_set(theme_bw())

## ---- opts ----
opts <- list()
opts$n <- 30 # samples
opts$p <- 250 # OTUs
opts$sigma <- .5 # noise
opts$k <- 6 # underlying OTU GPs
opts$sigma2_t <- .05
opts$Kf <- function(x, y) { (1 / (sqrt(2 * pi * opts$sigma2_t))) * exp(- 1 / (2 * opts$sigma2_t) * (x - y) ^ 2) }

## ---- make-sources ----
times <- runif(opts$n, 0, 1)

# underlying GPs

K <- matrix(0, opts$n, opts$n)
for(i in seq_len(opts$n)) {
  for(j in seq_len(opts$n)) {
    K[i, j] <- opts$Kf(times[i], times[j])
  }
}

Keig <- eigen(K)
Keig$values[Keig$values < 1e-10] <- 0
Ksqrt <- Keig$vectors %*% diag(sqrt(Keig$values)) %*% t(Keig$vectors)
Mu <- Ksqrt %*% matnorm(opts$n, opts$k, 1)
Mu_df <- data.frame(times = times, Mu)
mMu <- melt(Mu_df, id.vars = "times")

ggplot(mMu) +
  geom_point(aes(x = times, y = value, col = as.factor(variable))) +
  ggtitle("Latent Sources")

## ---- make-scores ----
# scores are clustered
memberships <- sample(opts$k, opts$p, replace = T)
X <- Mu[, memberships] + matnorm(opts$n, opts$p, opts$sigma)

mX <- melt(data.frame(times = times, X = X), id.vars = "times")
mX$variable <- as.numeric(str_extract(mX$variable, "[0-9]+"))
mX$membership <- as.factor(memberships[mX$variable])

ggplot(mX %>% filter(variable < 20)) +
  geom_point(aes(x = times, y = value, col = membership)) +
  facet_wrap(~variable)

## ---- pca ----
pca_res <- princomp(t(X))

V <- pca_res$loadings
class(V) <- "matrix"
V <- data.frame(times = times, V)
mV <- melt(V[, 1:10], id.vars = "times")
ggplot(mV) +
  geom_point(aes(x = times, y = value, col = variable)) +
  facet_wrap(~variable) +
  ggtitle("First 10 Loadings")

U <- pca_res$scores[, 1:10]
mU <- melt(data.frame(ix = 1:nrow(U), memberships = memberships, U),
           id.vars = c("memberships", "ix"))
mU$variable <- str_extract(mU$variable, "[0-9]+") %>%
  as.numeric()

mU_cast <- mU %>% dcast(ix + memberships ~ variable, value.var = "value")
colnames(mU_cast)[3:ncol(mU_cast)] <- paste0("X", colnames(mU_cast)[3:ncol(mU_cast)])

scores_plots <- list()
for(i in seq_len(5)) {
  for(j in seq_len(i - 1)) {
    scores_plots[[paste0(i, j)]] <- ggplot(mU_cast) +
      geom_point(aes_string(x = paste0("X", i), y = paste0("X", j), col = as.factor(memberships))) +
      ggtitle(sprintf("Scores, axes %d vs. %d", i, j))
  }
}

scores_plots
