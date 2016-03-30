
################################################################################
# A nonlinear multitable example, made from two concentric circles.
################################################################################

## ---- libraries ----
library("ggplot2")
theme_set(theme_bw())

## ---- generate-data ----
n <- 100 # number of samples
z <- as.factor(sample(0:1, n, replace = TRUE)) # class labels
X <- list() # will store combined data frames

r <- setNames(c(1, 1.5), c("0", "1")) # radii
X[[1]] <- cbind(r[z], 2 * pi * runif(n)) + matrix(runif(2 * n), n, 2)
X[[2]] <- cbind(r[z], 2 * pi * runif(n)) + matrix(runif(2 * n), n, 2)

#X[[2]] <- r[z] %*% t(c(1, 1)) + matrix(runif(2 * n), n, 2)

## ---- raw-data-plots ----
X1_df <- data.frame(X[[1]], z)
X1_df$x <- X1_df$X1 * cos(X1_df$X2)
X1_df$y <- X1_df$X1 * sin(X1_df$X2)

X2_df <- data.frame(X[[2]], z)
X2_df$x <- X2_df$X1 * cos(X2_df$X2)
X2_df$y <- X2_df$X1 * sin(X2_df$X2)

ggplot(X1_df) +
  geom_point(aes(x = x, y = y, col = z))
ggplot(X2_df) +
  geom_point(aes(x = x, y = y, col = z))

## ---- usual-cca ----
cancor_res <- cancor(X[[1]], X[[2]])
cancor_scores <- list(X1 = X[[1]] %*% cancor_res$xcoef,
                      X2 = X[[2]] %*% cancor_res$ycoef)

## ---- cca-plots ----
ggplot(data.frame(cancor_scores[[1]], z)) +
  geom_point(aes(x = X1, y = X2, col = z))
ggplot(data.frame(cancor_scores[[2]], z)) +
  geom_point(aes(x = X1, y = X2, col = z)) # don't even need second component

ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X1.1, col = z)) # perfect correlation between scores
ggplot(data.frame(cancor_scores, z)) + 
  geom_point(aes(x = X1.2, y = X1.2, col = z))
