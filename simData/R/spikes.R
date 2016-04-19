
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
  stopifnot(length(x0) == length(h) && length(x0) == length(w))
  function(x) {
    fx <- matrix(0, length(x), length(x0))
    for (j in seq_along(x0)) {
      fx[, j] <- symm_spike(x0[j], h[j], w[j])(x)
    }
    rowSums(fx)
  }
}
