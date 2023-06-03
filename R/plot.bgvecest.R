#' Plotting Draws of a VECX Submodel of a GVAR Model
#' 
#' A plot function for objects of class \code{"bgvecest"}.
#' 
#' @param x an object of class \code{"bgvecest"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param model indicates, for which model plots should be produced. Can be a character
#' vector of country names or an integer vector of the positions of the elements in
#' argument \code{x}.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.bgvecest <- function(x, ci = 0.95, type = "hist", model = NULL, ...) {
  
  if (is.null(model)) {
    
    for (i in 1:length(x)) {
      plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...)
    }
    
  } else {
    
    if (class(model) == "character") {
      model <- which(names(x) %in% model)
      if (length(model) == 0) {
        stop("There is no output for the specified country.")
      }
    }
    
    if (any(class(model) %in% c("numeric", "integer"))) {
      for (i in model) {
        plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...)
      }
    }
  }
}

