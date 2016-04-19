
################################################################################
# Functions for simulating spiky data
################################################################################

#' @title Functional for creating symmetric spikes
#' @description Returns functions that can be evaluated on x and returns f(x)
#' where f is zero outside of an interval, and is a symmetric spike within that
#' interval.
#' @param x0 The centerpoint for the peak.
#' @param h The maximum peak height.
#' @param w The half-width of the interval (the total width is 2 * w).
#' @examples
#' f <- symm_spike(10, 5, 3)
#' x <- seq(2, 20, .01)
#' plot(x, f(x), type = 'l')
#' f <- symm_spike(9, 3, 1)
#' lines(x, f(x))
#' f <- symm_spike(13, 1, 8)
#' lines(x, f(x))
#' @export
symm_spike <- function(x0, h, w) {
  function(x) {
    left_ix <- (x >= x0 - w) & (x <= x0)
    right_ix <- (x <= x0 + w) & (x >= x0)
    fx <- rep(0, length(x))
    fx[left_ix] <- (h / w) * (x[left_ix] - x0) + h
    fx[right_ix] <- - (h / w) * (x[right_ix] - x0) + h
    fx
  }
}

#' @title Functional for creating sum of symmetric spikes
#' @description Returns functions that can be evaluated on x and returns f(x)
#' where f is is a sum of symmetric spike functions.
#' @param x0 The centerpoints for the peak.
#' @param h The maximum peak heights.
#' @param w The half-widths of the intervals.
#' @export
#' @examples
#' x0 <- c(3, 15, 7)
#' h <- c(5, 3, 4)
#' w <- c(1, 4, .5)
#' f <- symm_spike_sum(x0, h, w)
#' x <- seq(2, 20, .01)
#' plot(x, f(x))
symm_spike_sum <- function(x0, h, w) {
  param_len <- c(length(x0), length(h), length(x0), length(w))
  stopifnot(abs(max(param_len) - min(param_len)) == 0)
  function(x) {
    fx <- matrix(0, length(x), length(x0))
    for (j in seq_along(x0)) {
      fx[, j] <- symm_spike(x0[j], h[j], w[j])(x)
    }
    rowSums(fx)
  }
}

#' @title Merge default parameters for random symmetric spike generation
#' @param opts A (partially filled) list of parameters to use when sampling
#' spikes. The currently available options are,
#'   $lambda_m The poisson parameter for the number of spikes to draw
#'   $lambda_h The exponential parameter for the height of the spikes.
#'   $lambda_w The exponential parameter for the width of the spikes.
#'   $alpha_x0 The first shape beta parameter for the positions of the
#' spikes.
#'   $beta_x0 The second shape beta parameter for the positions of the
#' spikes.
#'   $M The rescaling factor for the positions of the spikes (so they aren't
#'   just in the unit interval.
#' @export
merge_rsymm_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$lambda_m  <- 4
  default_opts$lambda_h <-  1
  default_opts$lambda_w  <- 1
  default_opts$alpha_x0  <- 1
  default_opts$beta_x0  <-  1
  default_opts$M <- 10
  modifyList(default_opts, opts)
}

#' @title Generate a list of random spike sum functions
#' @param n The number of random spike sum functions to return
#' @param opts A list specifying how the spike functions are created. See
#' merge_rsymm_opts() for available options.
#' @export
#' @examples
#' x <- seq(-4, 10, .01)
#' f_list <- rsymm_spike(10, list(M = 10)
#' fx <- sapply(f_list, function(f) { f(x) })
#' library("reshape2")
#' library("ggplot2")
#' mx <- melt(data.frame(x = x, fx = fx), id.vars = "x")
#' ggplot(mx) +
#'   geom_line(aes(x = x, y = value)) +
#'   facet_wrap(~variable) +
#'   theme_bw()
rsymm_spike <- function(n, opts = list()) {
  opts <- merge_rsymm_opts(opts)

  # generate spike parameters
  m <- rpois(n, opts$lambda_m)
  f_list <- list()
  for (i in seq_len(n)) {
    h <- rexp(m[i], opts$lambda_h)
    w <- rexp(m[i], opts$lambda_w)
    x0 <- opts$M * rbeta(m[i], opts$alpha_x0, opts$beta_x0)
    f_list[[i]] <- symm_spike_sum(x0, h, w)
  }
  f_list
}
