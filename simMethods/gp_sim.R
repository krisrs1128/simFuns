
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
