
################################################################################
# A nonlinear multitable example, made from two concentric circles.
################################################################################

## ---- kcca-libraries ----
# List of packages for session
.packages = c("kernlab",
              "plyr",
              "dplyr",
              "gridExtra",
              "ggplot2")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.packages[!.inst], repos='http://cran.rstudio.com/')
}

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
set.seed(04032016)

cat("\014")  # Clear console

rm(list=ls()) # Delete all existing variables
graphics.off() # Close all open plots

## ---- generate-data ----
n <- 100 # number of samples
z <- as.factor(sample(0:1, n, replace = TRUE)) # class labels
X <- list() # will store combined data frames

r <- setNames(c(1, 2.5), c("0", "1")) # radii
X[[1]] <- cbind(r[z], 2 * pi * rbeta(n, 1, 4)) + matrix(runif(2 * n), n, 2)
X[[2]] <- cbind(sqrt(r[z]), 2 * pi * runif(n)) + matrix(runif(2 * n), n, 2)
X <- X %>%
  lapply(function(x) {
    rownames(x) <- NULL
    x
  })

## ---- raw-data-plots ----
X1_df <- data.frame(X[[1]], z)
X1_df$x <- X1_df$X1 * cos(X1_df$X2)
X1_df$y <- X1_df$X1 * sin(X1_df$X2)

X2_df <- data.frame(X[[2]], z)
X2_df$x <- X2_df$X1 * cos(X2_df$X2)
X2_df$y <- X2_df$X1 * sin(X2_df$X2)

p1 <- ggplot(X1_df) +
  geom_point(aes(x = x, y = y, col = z))
p2 <- ggplot(X2_df) +
  geom_point(aes(x = x, y = y, col = z))
grid.arrange(p1, p2, ncol = 2)

## ---- usual-cca ----
cancor_res <- cancor(X[[1]], X[[2]])
cancor_scores <- list(X1 = X[[1]] %*% cancor_res$xcoef,
                      X2 = X[[2]] %*% cancor_res$ycoef)
cor(cancor_scores[[1]][, 1], cancor_scores[[2]][, 1]) %>%
  round(3)

## ---- cca-plots ----
p1 <- ggplot(data.frame(cancor_scores[[1]], z)) +
  geom_point(aes(x = X1, y = X2, col = z)) +
  ggtitle("CCA X-Scores (Concentric Circles)")
p2 <- ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z)) +
  ggtitle("CCA Correlation between X1, Y1 Scores")
grid.arrange(p1, p2, ncol = 2)

## ---- kcca ----
kcca_res <- kcca(X[[1]], X[[2]],
                 kernel = "rbfdot",
                 kpar = list(sigma = 2),
                 gamma = .25, ncomps = 2)
round(kcca_res@kcor, 3)
kcca_scores <- list(X1 = kcca_res@xcoef,
                    X2 = kcca_res@ycoef)

## ---- kcca-plots ----
p1 <- ggplot(data.frame(kcca_scores[[1]], z)) +
  geom_point(aes(x = X1, y = X2, col = z)) +
  ggtitle("KCCA X-Scores (Nonlinear Example)")
p2 <- ggplot(data.frame(kcca_scores)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z)) +
  ggtitle("KCCA Correlation between X1, Y1 Scores")
grid.arrange(p1, p2, ncol = 2)

## ---- fukumizu-example ----
# http://www.jmlr.org/papers/volume8/fukumizu07a/fukumizu07a.pdf
z <- runif(n)
U <- cbind(z %*% t(c(1 / sqrt(2), 1 / sqrt(2))) + .1 + matrix(rnorm(2 * n, 0, .05), n, 2),
           z %*% t(c(1 / sqrt(3), sqrt(2 / 3))) + 3 + matrix(rnorm(2 * n, 0, .05), n, 2))

R <- scale(U, center = F, scale = c(1.5, 1.5, 4.1, 4.1))
R <- sqrt(-5 * log(R))
pairs(R)

## ---- cca-fukumizu ----
cancor_res <- cancor(R[, 1:2], R[, 3:4])
cancor_scores <- list(X1 = R[, 1:2] %*% cancor_res$xcoef,
                      X2 = R[, 3:4] %*% cancor_res$ycoef)

## ---- cca-fukumizu-plots ----
ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X1.2, col = z))
ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X2.1, y = X2.2, col = z))

ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z)) # correlation between scores

## ---- kcca-fukumizu ----
kcca_res <- kcca(as.matrix(R[, 1:2]), as.matrix(R[, 3:4]),
                 kpar = list(sigma = 1))
kcca_scores <- list(X1 = kcca_res@xcoef, X2 = kcca_res@ycoef)

## ---- kcca-fukumizu-plots ----
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.1, y = X1.2, col = z))
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X2.1, y = X2.2, col = z))

ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z))
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.2, y = X2.2, col = z))

## ---- kcca-fukumizu2 ----
kcca_res <- kcca(as.matrix(R[, 1:2]), as.matrix(R[, 3:4]), kpar = list(sigma = .4))
kcca_scores <- list(X1 = kcca_res@xcoef, X2 = kcca_res@ycoef)

## ---- kcca-fukumizu-plots2 ----
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.1, y = X1.2, col = z))
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X2.1, y = X2.2, col = z))

# remember, don't expect linear correlation
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z))
ggplot(data.frame(kcca_scores, z)) +
  geom_point(aes(x = X1.2, y = X2.2, col = z))
