
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
