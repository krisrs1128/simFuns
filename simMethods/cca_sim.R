
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

## Standard Poisson formulation ------------------------------------------------

# Generate environmental data
X <- list()
X[["env"]] <- matnorm(opts$n, opts$p[2], opts$sigma_env)

# create species means and variances
mu_species <- rnorm(opts$p[1])
var_species <- rgamma(opts$p[1], 20)

plot(mu_species, ylim = c(-30, 30))
points(mu_species + 1.96 * var_species, col = 'red')
points(mu_species - 1.96 * var_species, col = 'red')

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

hist(exp(log_lambda_species[, 1]), 10)
hist(exp(log_lambda_species[, 2]), 10)
pairs(exp(log_lambda_species[, 1:5]))

X[["species"]] <- matrix(0, opts$n, opts$p[1])
for(k in seq_len(opts$p[1])) {
  X[["species"]][, k] <- rpois(opts$n, lambda = exp(log_lambda_species[, k]))
}

X[["species"]][rowSums(X[["species"]]) == 0, 1] <- 1

hist(X[["species"]][, 1])
hist(X[["species"]][, 2])
pairs(X[["species"]][, 1:10])
image(log(1 + X[["species"]]))
