
################################################################################
# Clustered regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] for r = 1:p[1] is a noisy realization of a much smaller number of
# b[1], ..., b[k]
################################################################################

## ---- libraries ----
library("simData")
library("expm")
library("plyr")
library("reshape2")
library("bayesMult")
library("ggplot2")
library("dplyr")
theme_set(theme_bw())

## ---- opts ----
opts <- list()
opts$n <- 200
opts$p <- c(100, 2) # 60 response dimensions, 20 features
opts$k <- 3
opts$sigma_b <- 1
opts$sigma_y <- 4

## ---- fixed-params ----
S <- matrix(c(0, 2, 3, 1, 2, 4), nrow = 2)
Psi <- diag(opts$p[2])

## ---- random-params ----
W <- sample(1:opts$k, opts$p[1], replace = T) # each equally likely
dfw <- data.frame(w = as.factor(W))
W <- t(model.matrix(~ w - 1, dfw, intercept = F))

E <- matnorm(opts$p[2], opts$p[1], opts$sigma_b)
B <- S %*% W + sqrtm(Psi) %*% E

## ---- generate-data ----
x_list <- list()
y_list <- list()

unique_ts <- list()
for(p in seq_len(opts$p[2])) {
  unique_ts[[p]] <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = opts$n)[-1]
}

for(r in seq_len(opts$p[1])) {
  x_list[[r]] <- matrix(0, opts$n, opts$p[2])
  for(p in seq_len(opts$p[2])) {
    x_list[[r]][, p] <- unique_ts[[p]]
  }
  y_list[[r]] <- x_list[[r]] %*% B[, r] + rnorm(opts$n, 0, opts$sigma_y)
}

## ---- plot-x ----
mX <- melt(x_list[[1]])
ggplot(mX) +
  geom_line(aes(x = Var1, y = value, col = as.factor(Var2))) +
  ggtitle("Two underlying OTUs")

## ---- plot-y ----
mY <- melt(y_list)
mB <- melt(B) %>%
  dcast(Var2 ~ Var1)
colnames(mB) <- c("L1", "beta1", "beta2")
mY <- mY %>%
  left_join(mB)
ggplot(mY) +
  geom_line(aes(x = Var1, y = value, group = L1, col = beta1)) +
  ggtitle("Genes, colored by coef. first OTU")
ggplot(mY) +
  geom_line(aes(x = Var1, y = value, group = L1, col = beta2)) +
  ggtitle("Genes, colored by coef. second OTU")

merged_lists <- lapply(seq_along(x_list), function(i) {
  data.frame(t = 1:opts$n, y = y_list[[i]], x = x_list[[i]])
})

m_merged <- melt(merged_lists, id.vars = "t") %>%
  filter(L1 < 10)

clust <- as.factor(apply(W, 2, which.max))
m_merged$cluster <- clust[m_merged$L1]
ggplot(m_merged) +
  geom_line(aes(x = t, y = value, linetype = variable, col = cluster)) +
  facet_wrap(~L1) +
  ggtitle("Clustered Series")

## ---- run-vb ----
phi <- rep(1 / opts$k, opts$k)
sigma <- opts$sigma_y # assume we know sigma for now

S0 <- matrix(c(1.5, 2.5, 2.5, 1.5, .5, 2.0), nrow = 2)
Psi0 <- diag(runif(opts$p[2], .8, 1.2))
data_list <- list(x_list = x_list, y_list = y_list)
param_list <- list(S = S0, Psi = Psi0, phi = phi, sigma = sigma)

vb_res <- vb_multinom(data_list, param_list, 4)

pi_guess <- t(vb_res$var_list$pi)
m <- t(vb_res$var_list$m)
clust <- as.factor(apply(pi_guess, 1, which.max))
m_df <- data.frame(m, clust = clust)
colnames(m_df) <- c("beta1", "beta2", "cluster")

## ---- plot-clusters ----
ggplot(m_df) +
  geom_point(aes(x = beta1, y = beta2, col = cluster)) +
  ggtitle("Recovered Clusters")

## ---- coord-desc ----
qplot(seq_along(vb_res$elbo), vb_res$elbo) +
  ggtitle("Objective after each update")

## ---- naive-model ----
B_lm <- matrix(0, nrow(B), ncol(B))
for(r in seq_len(ncol(B))) {
  B_lm[, r] <- coef(lm(y_list[[r]] ~ -1 + x_list[[r]]))
}

## ---- plot-b-hat ----
B_hat <- vb_res$param_list$S %*% vb_res$var_list$pi
mB_hat <- melt(B_hat)
colnames(mB_hat) <- c("var", "L1", "beta_hat")

mB_lm_hat <- melt(B_lm)
colnames(mB_lm_hat) <- c("var", "L1", "beta_lm_hat")

mB <- melt(B)
colnames(mB) <- c("var", "L1", "beta")
mB <- data.frame(mB, beta_lm_hat = mB_lm_hat$beta_lm_hat,
                 beta_hat = mB_hat$beta_hat)

merged_cast <- m_merged %>%
  dcast(t + L1 + cluster ~ variable) %>%
  left_join(mB) %>%
  filter(var == "1", L1 <= 9) %>%
  arrange(L1, t)

ggplot(merged_cast) +
  geom_point(aes(x = x.1, y = y, col = cluster), size = .5) +
  geom_abline(aes(intercept = 0, slope = beta), col = 'purple') +
  geom_abline(aes(intercept = 0, slope = beta_lm_hat), col = 'steelblue') +
  geom_abline(aes(intercept = 0, slope = beta_hat), col = 'indianred') +
  facet_wrap(~L1) +
  ggtitle("OTU 1 vs. 9 Genes")
