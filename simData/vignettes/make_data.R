
################################################################################
# Functions for making multitable sim data
################################################################################

## ---- setup ----
library("simData")
path <- "~/Documents/programming/simFuns/"
setwd(path)

## pls-data ----
pls_two_table <- ppls_shared(list(n = 5000, p = 40))
pairs(pls_two_table$X[[1]][, c(1:5)]) # first five columns of first table
pairs(cbind(pls_two_table$X[[1]][, c(1:4)], pls_two_table$X[[2]][, c(1:4)])) # first five columns, across first two tables
use_data(pls_two_table, pkg = "simData", overwrite = T)

pls_four_table <- ppls_shared(list(n = 5000, p = 40, k_unique = rep(2, 4)))
use_data(pls_four_table, pkg = "simData", overwrite = T)
