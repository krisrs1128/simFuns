
################################################################################
# Simulate canonical correspondence, according to models in
# https://web.stanford.edu/~hastie/Papers/zhulda.pdf
################################################################################

## ---- packages ----
library("simData")
library("vegan")

## ---- opts ----
opts <- list()
opts$n <- 40 # number of sites
opts$p <- c(100, 10) # number of species and environmental variables
opts$sigma_env <- 10

# Standard Poisson formulation ------------------------------------------------

## ---- generate-data ----
# Generate environmental data
X <- list()
X[["env"]] <- matnorm(opts$n, opts$p[2], opts$sigma_env)

# create species means and variances
mu_species <- rnorm(opts$p[1])
var_species <- rgamma(opts$p[1], 20)

# define a new gradient
alpha <- rnorm(opts$p[2])
alpha <- alpha / sqrt(sum(alpha ^ 2))

# gaussian response curves baseline
a <- rnorm(opts$p[1], log(sqrt(2 * pi)), 1)

# generate species counts
log_lambda_species <- matrix(0, opts$n, opts$p[1])

for(i in seq_len(opts$n)) {
  for(k in seq_len(opts$p[1])) {
    log_lambda_species[i, k] <- a[k] - (X[["env"]][i, ] %*% alpha - mu_species[k]) ^ 2 / var_species[k]
  }
}

X[["species"]] <- matrix(0, opts$n, opts$p[1])
for(k in seq_len(opts$p[1])) {
  X[["species"]][, k] <- rpois(opts$n, lambda = exp(log_lambda_species[, k]))
}

X[["species"]][rowSums(X[["species"]]) == 0, 1] <- 1

## ---- plot-data-sim ----
plot(mu_species, ylim = c(-50, 50))
points(mu_species + 1.96 * var_species, col = 'red')
points(mu_species - 1.96 * var_species, col = 'red')

hist(exp(log_lambda_species[, 1]), 10)
hist(exp(log_lambda_species[, 2]), 10)
pairs(exp(log_lambda_species[, 1:5]))

hist(X[["species"]][, 1])
hist(X[["species"]][, 2])
pairs(X[["species"]][, 1:10])
pairs(log(1 + X[["species"]][, 1:10]))
image(log(1 + X[["species"]]))

## ---- apply-cca ----
cca_res <- cca(X[[2]], X[[1]])
plot(cca_res)

## ---- evaluate-scores ----
site_scores <- cca_res$CCA$u[, 1:2]
env_proj <- X[["env"]] %*% alpha

ggplot(data.frame(site_scores, env = env_proj)) +
  geom_point(aes(x = CCA1, y = CCA2, col = env)) +
  scale_color_gradient2(mid = "antiquewhite") +
  ggtitle("Recovered Environment Scores vs. Truth")

species_scores <- cca_res$CCA$v[, 1:2]
ggplot(data.frame(species_scores, mu = mu_species, sigma2 = var_species)) +
  geom_point(aes(x = CCA1, y = CCA2, col = mu, size = sigma2)) +
  scale_color_gradient2(mid = "antiquewhite") +
  ggtitle("Recovered Species Scores vs. Species Means")

ggplot(data.frame(species_scores, a = a, sigma2 = var_species)) +
  geom_point(aes(x = CCA1, y = CCA2, col = a, size = sigma2)) +
  scale_color_gradient2(mid = "antiquewhite") +
  ggtitle("Recovered Species Scores vs. Species Baseline Abundance")

## ---- cancor ----
cancor_res <- CCorA(log(1 + X[["species"]]), X[["env"]])

cancor_env <- list()
cancor_env[[1]] <- data.frame(cancor_res$Cy, type = "Cy",
                         env_proj = X[["env"]] %*% alpha)
cancor_env[[2]] <- data.frame(cancor_res$Cx, type = "Cx",
                         env_proj = X[["env"]] %*% alpha)
cancor_env <- do.call(rbind, cancor_env)

ggplot(cancor_env) +
  geom_point(aes(x = CanAxis1, y = CanAxis2, col = env_proj)) +
  facet_wrap(~type)

cancor_species <- list()
cancor_species[[1]] <- data.frame(cancor_res$corr.Y.Cy, type1 = "Y",
                                  type2 = "Cy", species_means = mu_species,
                                  species_vars = var_species)
cancor_species[[2]] <- data.frame(cancor_res$corr.Y.Cx, type1 = "Y",
                                  type2 = "Cx", species_means = mu_species,
                                  species_vars = var_species)
cancor_species <- do.call(rbind, cancor_species)
ggplot(cancor_species) +
  geom_point(aes(x = CanAxis1, y = CanAxis2, col = species_means, size = species_vars)) +
  facet_grid(type1 ~ type2)
