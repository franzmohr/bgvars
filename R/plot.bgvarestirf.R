#' Plotting Impulse Responses of Bayesian GVAR Submodels
#' 
#' A plot function for objects of class \code{"bgvarestirf"}.
#' 
#' @param x an object of class \code{"bgvarestirf"}, usually, a result of a call to \code{\link{irf.bgvarest}}.
#' @param ... further graphical parameters.
#' 
#' @export 
#' @rdname irf.bgvarest
plot.bgvarestirf <- function(x, ...) {
  for (i in 1:length(x)) {
    plot(x[[i]], ...)
  }
}


