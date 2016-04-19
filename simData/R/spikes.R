
################################################################################
# Functions for simulating spiky data
################################################################################

#' @title Functional for creating spikes
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
