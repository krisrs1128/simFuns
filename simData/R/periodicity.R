
################################################################################
# Functions to generate data with periodicities
################################################################################

#' @title Functional returning a periodic function
#' @description Generates the function
#' h * sin(2 * pi * x * period) I(x \in [x_left x_right])
#' @param x_left Left endpoint of the nonzero interval.
#' @param x_right Right endpoint of the nonzero interval.
#' @param h The maximum height of the oscillation.
#' @param period The period of the function.
#' @examples
#' f1 <- periodicity(4, 8, 3, 2)
#' x <- seq(2, 10, .01)
#' plot(x, f1(x), type = 'l')
#' f2 <- periodicity(5, 9, 1, .6)
#' lines(x, f2(x))
#' plot(x, cumsum(2 + f1(x) + f2(x)), type = "l")
#' @export
periodicity <- function(x_left, x_right, h, period) {
  function(x) {
    interval_ix <- (x >= x_left) & (x <= x_right)
    fx <- rep(0, length(x))
    fx[interval_ix] <- h * sin(x[interval_ix] * (2 * pi * period))
    fx
  }
}

#' @title Functional returning sum of periodic functions
#' @description Generates the function
#' \sum_{j} h[j] * sin(2 * pi * x * period[j]) I(x \in [x_left[j] x_right[j]])
#' @param x_left Left endpoints of the nonzero intervals.
#' @param x_right Right endpoints of the nonzero intervals.
#' @param h The maximum heights of the oscillations.
#' @param period The periods of the functions.
#' @examples
#' x <- seq(2, 10, .01)
#' x_left <- c(1, 5, 6)
#' x_right <- c(3, 9, 7)
#' h <- c(2, 1, 3)
#' period <- c(1, 5, 2)
#' f <- periodicity_sum(x_left, x_right, h, period)
#' plot(x, f(x), type = "l")
#' plot(x, cumsum(1 + f(x)), type = "l")
#'@export
periodicity_sum <- function(x_left, x_right, h, period) {
  param_len <- c(length(x_left), length(x_right), length(h), length(period))
  stopifnot(abs(max(param_len) - min(param_len)) == 0)
  function(x) {
    fx <- matrix(0, length(x), length(x_left))
    for (j in seq_along(x_left)) {
      fx[, j] <- periodicity(x_left[j], x_right[j], h[j], period[j])(x)
    }
    rowSums(fx)
  }
}

#' @title Merge default options when generating periodic functions
merge_rperiodicity_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$lambda_m  <- 4
  default_opts$lambda_h <- 1
  default_opts$alpha_per <- 2
  default_opts$lambda_per <- 3
  default_opts$lambda_w  <- 1
  default_opts$alpha_x0  <- 1
  default_opts$beta_x0  <-  1
  default_opts$M <- 10
  modifyList(default_opts, opts)
}

#' @title Generate functions that are random sums of periodicities 
#' @param n The total number of functions to return.
#' @param opts A list of options for how to generate the periodicity parameters.
#' See merge_rperiodicity_opts() for possible arguments.
#' @return A list of functions with random parameters, which can be evaluated
#' on new x's.
#' @examples
#' x <- seq(0, 10, .01)
#' f_list <- rperiodicity(10, list(M = 10))
#' fx <- sapply(f_list, function(f) { f(x) })
#' library("reshape2")
#' library("ggplot2")
#' mx <- melt(data.frame(x = x, fx = fx), id.vars = "x")
#' ggplot(mx) +
#'   geom_line(aes(x = x, y = value)) +
#'   facet_wrap(~variable) +
#'    theme_bw()
#' @export
rperiodicity <- function(n, opts = list()) {
  opts <- merge_rperiodicity_opts(opts)

  # generate spike parameters
  m <- rpois(n, opts$lambda_m)
  f_list <- list()
  for (i in seq_len(n)) {
    h <- rexp(m[i], opts$lambda_h)
    w <- rexp(m[i], opts$lambda_w)
    x0 <- opts$M * rbeta(m[i], opts$alpha_x0, opts$beta_x0)
    period <- rgamma(m[i], opts$alpha_per, opts$lambda_per)
    f_list[[i]] <- periodicity_sum(x0, x0 + w, h, period)
  }
  f_list
}
