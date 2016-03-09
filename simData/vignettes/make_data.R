
################################################################################
# Functions for making multitable sim data
################################################################################

## setwd ----
path <- "~/Documents/programming/simFuns/"
setwd(path)

## pls-data ----
pls_two_table <- ppls_shared(list(n = 5000, p = 40))
pairs(pls_data$X[[1]][, c(1:5)]) # first five columns of first table
pairs(cbind(pls_data$X[[1]][, c(1:4)], pls_data$X[[2]][, c(1:4)])) # first five columns, across first two tables

pls_four_table <- ppls_shared(list(n = 5000, p = 40, k_unique = rep(2, 4)))
plot(pls_data$X[[1]][, 1:2]) # first two columns of first table
use_data(pls_data, pkg = "simData")
