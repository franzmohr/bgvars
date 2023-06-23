#' Plotting Draws of a Bayesian Global VAR Models
#' 
#' A plot function for objects of class \code{"bgvarest"}.
#' 
#' @param x an object of class \code{"bgvarestirf"}, usually, a result of a call to \code{\link{irf.bgvarest}}.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.bgvarestirf <- function(x, ...) {
  for (i in 1:length(x)) {
    plot(x[[i]], ...)
  }
}


