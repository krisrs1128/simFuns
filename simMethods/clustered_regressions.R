
################################################################################
# Clustered regression models, a special case of
# "Flexible latent variable models for multi-task learning"
#
# b[r] for r = 1:p[1] is a noisy realization of a much smaller number of
# b[1], ..., b[k]
################################################################################

## ---- cr-setup ----
# List of packages for session
.packages = c("expm",
              "plyr",
              "reshape2",
              "bayesMult",
              "ggplot2",
              "dplyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.packages[!.inst], repos='http://cran.rstudio.com/')
}

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
set.seed(5052016)

cat("\014")  # Clear console

rm(list=ls()) # Delete all existing variables
graphics.off() # Close all open plots

theme_set(theme_bw())
min_theme <- theme_update(panel.border = element_blank(),
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

## ---- fixed-params ----
S <- matrix(rnorm(opts$p[2] * opts$k, 0, sd = opts$sigma_s),
            opts$p[2], opts$k)
Psi <- diag(opts$p[2])

## ---- random-params ----
W <- sample(1:opts$k, opts$p[1], replace = T) # each equally likely
dfw <- data.frame(w = as.factor(W))
W <- t(model.matrix(~ w - 1, dfw, intercept = F))
true_clust <- setNames(apply(W, 2, which.max), 1:opts$p[1])
W[, 1:10]

E <- matrix(rnorm(prod(opts$p), sd = opts$sigma_b), opts$p[2], opts$p[1])
B <- S %*% W + sqrtm(Psi) %*% E

## ---- plot-B ----
m_b <- melt(B)
m_b$true_clust <- true_clust[m_b$Var2]
m_b$Var1 <- factor(m_b$Var1, opts$p[1]:1)

ggplot(m_b) +
  geom_tile(aes(y = Var1, x = reorder(Var2, true_clust), fill = value)) +
  scale_fill_gradient2(midpoint = 0, high = "#90ee90", low = "#000080") +
  labs(fill = "Coef. Value", x = "Tasks [reordered]", y = "Features") +
  theme(axis.text.x = element_text(angle = -90))

## ---- generate-data ----
x_list <- list()
y_list <- list()

for(r in seq_len(opts$p[1])) {
  x_list[[r]] <- matrix(rnorm(opts$n * opts$p[2]), opts$n, opts$p[2])
  y_list[[r]] <- x_list[[r]] %*% B[, r] + rnorm(opts$n, 0, opts$sigma_y)
}

## ---- join-data ----
m_y <- melt(y_list, value.name = "y")
m_y$Var2 <- NULL
colnames(m_y) <- c("time", "y", "task")
m_x <- melt(x_list, value.name = "x")
colnames(m_x) <- c("time", "feature", "x", "task")
colnames(m_b) <- c("feature", "task", "slope", "cluster")
m_b$feature <- m_b$feature %>%
  as.character %>%
  as.numeric

merged_data <- m_x %>%
  left_join(m_y) %>%
  left_join(m_b)

## ---- vis-reg-lines ----
ggplot(merged_data %>%
         filter(task < 11, feature < 6)) +
  geom_point(aes(x = x, y = y, col = as.factor(cluster)), size = .3) +
  geom_abline(aes(slope = slope, intercept = 0, col = as.factor(cluster))) +
  scale_x_continuous(breaks = c(-10, 0, 10)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  labs(col = "Cluster") +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67")) +
  facet_grid(feature ~ task)

## ---- run-vb ----
phi <- rep(1 / opts$k, opts$k)
sigma <- opts$sigma_y # assume we know sigma for now

S0 <- matrix(rnorm(opts$p[2] * opts$k, 0, sd = opts$sigma_s),
            opts$p[2], opts$k)
Psi0 <- diag(runif(opts$p[2], .5, 1.5))
data_list <- list(x_list = x_list, y_list = y_list)
param_list <- list(S = S0, Psi = Psi0, phi = phi, sigma = sigma)
vb_res <- vb_multinom(data_list, param_list, 1)
str(vb_res)

## ---- plot-pi ----
pi_guess <- t(vb_res$var_list$pi)
true_clust <- apply(W, 2, function(x) which(x == 1))
pi_est <- melt(data.frame(pi = pi_guess, true_clust),
               measure.vars = paste("pi", 1:opts$k, sep = "."))
pi_est$true  <- as.factor(pi_est$true)

ggplot(pi_est) +
  geom_jitter(aes(x = variable, y = value, col = true)) +
  labs(col = "True Cluster", x = "Component", y = "Estimated Probability") +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67"))

## ---- table-clusters ----
clust <- as.factor(apply(pi_guess, 1, which.max))
table(clust, true_clust)

## ---- naive-model ----
B_lm <- matrix(0, nrow(B), ncol(B))
V <- array(0, c(nrow(B), nrow(B), ncol(B)))
for(r in seq_len(ncol(B))) {
  lm_fit <- lm(y_list[[r]] ~ -1 + x_list[[r]])
  B_lm[, r] <- coef(lm_fit)
  V[,, r] <- vcov(lm_fit)
}

## ---- merge-fits ----
m_vb <- melt(t(vb_res$var_list$m))
m_lm <- melt(t(B_lm))
m_vjj_lm <- apply(V, 3, diag) %>%
  t() %>%
  melt()
m_vjj_vb <- apply(vb_res$var_list$v, 3, diag) %>%
  t() %>%
  melt()

colnames(m_vb) <- c("task", "feature", "vb_slope")
colnames(m_lm) <- c("task", "feature", "lm_slope")
colnames(m_vjj_lm) <- c("task", "feature", "lm_var")
colnames(m_vjj_vb) <- c("task", "feature", "vb_var")

merged_data_fit  <- merged_data %>%
  left_join(m_vb) %>%
  left_join(m_lm) %>%
  left_join(m_vjj_lm) %>%
  left_join(m_vjj_vb)

m_merged_data_fit <- merged_data_fit %>%
  melt(measure.vars = c("vb_slope", "lm_slope"),
       variable.name = c("slope_type"), value.name = "slope_est") %>%
  melt(measure.vars = c("lm_var", "vb_var"), variable.name = "var_type",
       value.name = c("variance")) %>%
  filter(!(slope_type == "lm_slope" & var_type == "vb_var"),
         !(slope_type == "vb_slope" & var_type == "lm_var"))

## ---- plot-betas ----
plot_betas_df <- m_merged_data_fit %>%
  filter(time == 1)

methods_label <- c("vb_slope" = "Bayesian Multitask",
                   "lm_slope" = "Independent Regressions")
ggplot(plot_betas_df) +
  geom_abline(slope = 1, intercept = 0, col = "#db4551", size = .5) +
  geom_errorbar(aes(x = slope, ymin = slope_est - 1.9 * sqrt(variance),
                    ymax = slope_est + 1.9 * sqrt(variance)),
                size = .05, col = "#292929") +
  geom_point(aes(x = slope, y = slope_est), col =  "#292929", size = .2) +
  coord_fixed() +
  labs(color = "Method", y = "Estimated Slope", x = "True Slope") +
  facet_grid(. ~ slope_type, labeller = as_labeller(methods_label))

## ---- plot-fitted-reg ----
methods_label_short <- c(c("vb_slope" = "BM",
                           "lm_slope" = "LM"),
                         setNames(1:10, 1:10))
ggplot(m_merged_data_fit %>%
         filter(task < 10, feature < 5)) +
  geom_point(aes(x = x, y = y, col = as.factor(cluster)), size = .3) +
  geom_abline(aes(slope = slope_est, intercept = 0, linetype = slope_type,
                  col = as.factor(cluster))) +
  geom_abline(aes(slope = slope_est + 1.9 * sqrt(variance), intercept = 0,
                  col = as.factor(cluster), linetype = slope_type), alpha = 0.4) +
  geom_abline(aes(slope = slope_est - 1.9 * sqrt(variance), intercept = 0,
                  col = as.factor(cluster), linetype = slope_type), alpha = 0.4) +
  geom_abline(aes(slope = slope, intercept = 0, linetype = "True",
                  col = as.factor(cluster))) +
  scale_y_continuous(breaks = c(-15, 0, 15)) +
  scale_color_manual(values = c("#5FABC8", "#ffdead", "#c16a67")) +
  scale_linetype(labels = c("True", "LM", "BM")) +
  labs(linetype = "Slope", col = "Task Cluster") +
  facet_grid(feature ~ slope_type ~ task,
             labeller = as_labeller(methods_label_short))
