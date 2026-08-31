#' Plotting Draws of Submodels of a GVAR Model
#' 
#' A plot function for objects of class 'submodelestlist'.
#' 
#' @param x an object of class 'submodelestlist', usually, a result of a call to
#' \code{\link[bvartools]{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param variables character vector of variables that should be plotted. Default is \code{"all"}.
#' Other options are \code{"domestic"}, \code{"foreign"}, \code{"global"}, \code{"deterministic"}
#' and \code{"sigma"}.
#' @param group character. Name of the element in argument \code{x}, for which posterior draws
#' should be plotted. If \code{NULL} (default), all submodels are used.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.submodelestlist <- function(x, ci = 0.95, type = "hist", variables = "all", group = NULL, ...) {
  
  pos <- 1:length(x)
  if (!is.null(group)) {
    pos <- which(names(x) %in% group)
    if (length(pos) == 0) {
      stop("There is no output for the specified group.")
    } 
  }
  
  for (i in pos) {
    plot(x[[i]], ci = ci, type = type, variables = variables, ctry = names(x)[i], ...)
  }
  
}


