#' Plotting Draws of a Bayesian Global VAR Models
#' 
#' A plot function for objects of class \code{"bgvarest"}.
#' 
#' @param x an object of class \code{"bgvarest"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param ctry character. Name of the element in argument \code{x}, for which posterior draws
#' should be plotted. If \code{NULL} (default), all submodels are used.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.bgvarest <- function(x, ci = 0.95, type = "hist", ctry = NULL, ...) {
  
  pos <- 1:length(x)
  if (!is.null(ctry)) {
    pos <- which(names(x) %in% ctry)
    if (length(pos) == 0) {
      stop("There is no output for the specified country.")
    } 
  }
  
  for (i in pos) {
    plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...)
  }
  
}


