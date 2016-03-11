
################################################################################
# Multiresolution GP simulation
################################################################################


## ---- packages ----
library("ggplot2")
library("simData")
library("plyr")
library("dplyr")
library("reshape2")
library("stringr")
theme_set(theme_bw())

## ---- utils ----
Kf_gen <- function(sigma) {
  function(x, y) {
    (1 / (sqrt(2 * pi * sigma))) * exp(-1 / (2 * sigma) * (x - y) ^ 2)
  }
}

gen_K_subset <- function(times, lower, upper, sigma) {
  cur_ix <- which(times > lower & times < upper)
  cur_n <- length(cur_ix)
  K_sub <- matrix(0, cur_n, cur_n)
  for(i in seq_len(cur_n)) {
    for(j in seq_len(cur_n)) {
      K_sub[i, j] <- Kf_gen(sigma)(times[cur_ix[i]], times[cur_ix[j]])
    }
  }
  K_sub
}

## ---- opts ----
set.seed(2016)
opts <- list()
opts$n <- 70 # samples
opts$p <- 250 # OTUs
opts$sigma <- .5 # noise
opts$k <- 6 # underlying OTU GPs

## ---- make-sources ----
times <- runif(opts$n, 0, 1)

# underlying GPs
K <- array(0, dim = c(opts$n, opts$n, 2))
K[,, 1] <- gen_K_subset(times, 0, 1, 0.06)
K[times < .45, times < .45, 2] <- gen_K_subset(times, 0, .45, .03)
K[times > .55, times > .55, 2] <- gen_K_subset(times, .55, 1, .03)
K[times >= .45 & times <= .55, times >= .45 & times <= .55, 2] <- gen_K_subset(times, .45, .55, .02)

K <- apply(K, 1:2, sum) # sum over different resolutions
Keig <- eigen(K)
Keig$values[Keig$values < 1e-10] <- 0
Ksqrt <- Keig$vectors %*% diag(sqrt(Keig$values)) %*% t(Keig$vectors)
Mu <- Ksqrt %*% matnorm(opts$n, opts$k, 1)
Mu_df <- data.frame(times = times, Mu)
mMu <- melt(Mu_df, id.vars = "times")

ggplot(mMu) +
  geom_point(aes(x = times, y = value, col = as.factor(variable))) +
  ggtitle("MRA Latent Sources") +
  facet_wrap(~variable)

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
      ggtitle(sprintf("Scores, axes %d vs. %d", i, j)) +
      theme(legend.title = element_blank())
  }
}

scores_plots
