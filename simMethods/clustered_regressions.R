
################################################################################
# Clustered regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] for r = 1:p[1] is a noisy realization of a much smaller number of
# b[1], ..., b[k]
################################################################################

## ---- libraries ----
library("expm")
library("plyr")
library("reshape2")
library("bayesMult")
library("ggplot2")
library("dplyr")

theme_set(theme_bw())
min_theme <- theme(panel.border = element_blank(),
                   panel.grid = element_blank(),
                   axis.ticks = element_blank(),
                   legend.title = element_text(size = 10),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 10),
                   strip.background = element_blank(),
                   strip.text = element_text(size = 10),
                   legend.key = element_blank())

## ---- opts ----
opts <- list()
opts$n <- 15
opts$p <- c(100, 10)
opts$k <- 3
opts$sigma_b <- 1
opts$sigma_s <- 2
opts$sigma_y <- 2
set.seed(5052016)

## ---- fixed-params ----
S <- matrix(rnorm(opts$p[2] * opts$k, 0, sd = opts$sigma_s),
            opts$p[2], opts$k)
#S[1:(opts$p[2] / 2), ] <- S[1:(opts$p[2] / 2)] + 1
#S[(opts$p[2] / 2 + 1):opts$p[2], ] <- S[(opts$p[2] / 2 + 1):opts$p[2], ] - 1
#S <- S[sample(1:opts$p[2]), ]
Psi <- diag(opts$p[2])

## ---- random-params ----
W <- sample(1:opts$k, opts$p[1], replace = T) # each equally likely
dfw <- data.frame(w = as.factor(W))
W <- t(model.matrix(~ w - 1, dfw, intercept = F))

E <- matrix(rnorm(prod(opts$p), sd = opts$sigma_b), opts$p[2], opts$p[1])
B <- S %*% W + sqrtm(Psi) %*% E
true_clust <- setNames(apply(W, 2, which.max), 1:opts$p[1])

## ---- plot-B ----
m_b <- melt(B)
m_b$true_clust <- true_clust[m_b$Var2]

ggplot(m_b) +
  geom_tile(aes(y = Var1, x = reorder(Var2, true_clust), fill = value)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  min_theme +
  xlab("Tasks [reordered]") +
  ylab("Features") +
  theme(axis.text.x = element_text(angle = -90))

## ---- generate-data ----
x_list <- list()
y_list <- list()

unique_ts <- list()
for(p in seq_len(opts$p[2])) {
  #unique_ts[[p]] <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = opts$n)[-1]
  unique_ts[[p]] <- rnorm(opts$n)
}

for(r in seq_len(opts$p[1])) {
  x_list[[r]] <- matrix(0, opts$n, opts$p[2])
  for(p in seq_len(opts$p[2])) {
    x_list[[r]][, p] <- unique_ts[[p]]
  }
  y_list[[r]] <- x_list[[r]] %*% B[, r] + rnorm(opts$n, 0, opts$sigma_y)
}

## ---- join-data ----
m_y <- melt(y_list, value.name = "y")
m_y$Var2 <- NULL
colnames(m_y) <- c("time", "y", "task")
m_x <- melt(x_list, value.name = "x")
colnames(m_x) <- c("time", "feature", "x", "task")
colnames(m_b) <- c("feature", "task", "slope", "cluster")

merged_data <- m_x %>%
  left_join(m_y) %>%
  left_join(m_b)

## ---- vis-reg-lines ----
ggplot(merged_data %>%
         filter(task < 10, feature < 5)) +
  geom_point(aes(x = x, y = y, col = as.factor(cluster)), size = .3) +
  geom_abline(aes(slope = slope, intercept = 0, col = as.factor(cluster))) +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67")) +
  facet_grid(feature ~ task) +
  min_theme

## ---- run-vb ----
phi <- rep(1 / opts$k, opts$k)
sigma <- opts$sigma_y # assume we know sigma for now

S0 <- matrix(rnorm(opts$p[2] * opts$k, 0, sd = opts$sigma_s),
            opts$p[2], opts$k)
Psi0 <- diag(runif(opts$p[2], .5, 1.5))
data_list <- list(x_list = x_list, y_list = y_list)
param_list <- list(S = S0, Psi = Psi0, phi = phi, sigma = sigma)
vb_res <- vb_multinom(data_list, param_list, 1)

## ---- plot-pi ----
pi_guess <- t(vb_res$var_list$pi)
true_clust <- apply(W, 2, function(x) which(x == 1))
pi_est <- melt(data.frame(pi = pi_guess, true_clust),
               measure.vars = paste("pi", 1:opts$k, sep = "."))
pi_est$true  <- as.factor(pi_est$true)

ggplot(pi_est) +
  geom_jitter(aes(x = variable, y = value, col = true)) +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67")) +
  min_theme

## ---- plot-clusters ----
m <- t(vb_res$var_list$m)
clust <- as.factor(apply(pi_guess, 1, which.max))
m_df <- data.frame(m, clust = clust)
colnames(m_df) <- c("beta_1", "beta_2", "cluster")
table(clust, true_clust)

## ---- naive-model ----
B_lm <- matrix(0, nrow(B), ncol(B))
for(r in seq_len(ncol(B))) {
  B_lm[, r] <- coef(lm(y_list[[r]] ~ -1 + x_list[[r]]))
}

## ---- merge-fits ----
m_vb <- melt(t(vb_res$var_list$m))
m_lm <- melt(t(B_lm))
colnames(m_vb) <- c("task", "feature", "vb_slope")
colnames(m_lm) <- c("task", "feature", "lm_slope")
merged_data_fit  <- merged_data %>%
  left_join(m_vb) %>%
  left_join(m_lm)

m_merged_data_fit <- merged_data_fit %>%
  melt(measure.vars = c("slope", "vb_slope", "lm_slope"),
       variable.name = c("slope_type"), value.name = "slope")

## ---- plot-betas ----
plot_betas_df <- merged_data_fit %>%
  filter(time == 1) %>%
  melt(measure.vars = c("vb_slope", "lm_slope"))

ggplot(plot_betas_df) +
  geom_abline(slope = 1, intercept = 0, col = "#696969") +
  geom_point(aes(x = slope, y = value, col = variable), alpha = 0.5) +
  facet_wrap(~variable) +
  coord_fixed() +
  min_theme

## ---- plot-fitted-reg ----
ggplot(m_merged_data_fit %>%
         filter(task < 10, feature < 5)) +
  geom_point(aes(x = x, y = y, col = as.factor(cluster)), size = .3) +
  geom_abline(aes(slope = slope, intercept = 0, linetype = slope_type,
                  col = as.factor(cluster))) +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67")) +
  facet_grid(feature ~ task) +
  min_theme
