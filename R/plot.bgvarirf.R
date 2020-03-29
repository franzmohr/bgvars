#' Plotting Impulse Responses of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarirf".
#' 
#' @param x an object of class "bvarirf", usually, a result of a call to \code{\link{irf}}.
#' @param ... further graphical parameters.
#' 
#' @export
#' @rdname irf
plot.bgvarirf <- function(x, ...) {
  impulse <- attr(x, "impulse")
  response <- attr(x, "response")
  x <- cbind(0, x)
  ylab <- paste(response[1], ": ", response[2], sep = "")
  stats::plot.ts(x, plot.type = "single", lty = c(1, 2, 1, 2), ylab = ylab, ...)
}