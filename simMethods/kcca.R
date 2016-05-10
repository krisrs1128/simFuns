
################################################################################
# A nonlinear multitable example, made from two concentric circles.
################################################################################

## ---- kcca-libraries ----
# List of packages for session
.packages <- c("kernlab", "plyr", "dplyr", "scales", "grid",
               "gridExtra", "ggplot2", "GGally")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Load packages into session 
sapply(.packages, require, character.only = TRUE)
set.seed(04032016)
cat("\014")  # Clear console
rm(list = ls()) # Delete all existing variables
graphics.off() # Close all open plots

## ---- fukumizu-example ----
# http://www.jmlr.org/papers/volume8/fukumizu07a/fukumizu07a.pdf
m <- 1000
n <- 100
z <- runif(m)
u <- z + 0.06 + rnorm(m, 0, .035)
v <- z + 3 + rnorm(m, 0, .035)
u <- u[u > 0 & u < 1.5][1:n]
v <- v[v > 0 & v < 4.1][1:n]
R <- cbind(sqrt(-4 * log(u / 1.5)), sqrt(-4 * log(v / 4.1)))

theta <- list("x" = runif(n, 0, 2 * pi),
              "y" = runif(n, 0, 2 * pi))
X <- data.frame("x" = c(R[, 1] * cos(theta[["x"]]),
                        R[, 2] * cos(theta[["y"]])),
                "y" = c(R[, 1] * sin(theta[["x"]]),
                        R[, 2] * sin(theta[["y"]])),
                "cluster" = c(rep("0", n), rep("1", n)))

ggplot(X) +
  geom_point(aes(x = x, y = y, col = cluster), size = .5, alpha = 0.7) +
  scale_color_manual(values = c("#008080", "#A9A3D7"))

## ---- cca-fukumizu ----
x_list <- list(X = as.matrix(X[X$cluster == "0", 1:2]),
               Y = as.matrix(X[X$cluster == "1", 1:2]))
cancor_res <- cancor(x_list[[1]], x_list[[2]])
cancor_scores <- list(X1 = x_list[[1]] %*% cancor_res$xcoef,
                      X2 = x_list[[2]] %*% cancor_res$ycoef)

## ---- cca-fukumizu-plots ----
ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X1.2, col = z)) +
  coord_fixed(ratio = cancor_res$cor[2] / cancor_res$cor[1])
ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X2.1, y = X2.2, col = z)) +
  coord_fixed(ratio = cancor_res$cor[2] / cancor_res$cor[1])

 # correlation between scores
ggplot(data.frame(cancor_scores, z)) +
  geom_point(aes(x = X1.1, y = X2.1, col = z)) +
  coord_fixed()

## ---- kcca-fukumizu ----
kcca_res <- kcca(x_list[[1]], x_list[[2]], kpar = list(sigma = .01))
kcca_scores <- list(X1 = kcca_res@xcoef, X2 = kcca_res@ycoef)
plot(x_list[[1]][, 1], kcca_scores$X1[, 1])

## ---- kcca-fukumizu-plots ----
x_res <- data.frame(X, score = c(kcca_scores$X1[, 1], kcca_scores$X2[, 1]))

ggplot(x_res) +
  geom_point(aes(x = x, y = y, col = score, shape = cluster)) +
  scale_color_gradient(low = "#36191E", high = "#E3809C")

ggplot(data.frame(kcca_scores, R = R[, 1])) +
  geom_point(aes(x = X1.1, y = X2.1, col = R))
ggplot(data.frame(kcca_scores, R = R[, 1])) +
  geom_point(aes(x = X1.2, y = X2.2, col = R))
