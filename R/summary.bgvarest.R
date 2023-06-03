#' Summarising Country-Specific VARX Submodels of a GVAR Model
#'
#' summary method for class \code{"bgvarest"}.
#'
#' @param object an object of class \code{"bgvarest"}, usually, a result of a call to
#' \code{\link{draw_posterior}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param period integer. Index of the period, for which the summary statistics should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param ctry character. Name of the element in argument \code{object}, for which summary
#' statistics should be obtained. If \code{NULL} (default), all country submodels are used.
#' @param x an object of class \code{"summary.bgvarest"}, usually, a result of a call to
#' \code{\link{summary.bgvarest}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.bgvarest} returns a list of class \code{"summary.bgvarest"},
#' which contains summary statistics of country-specific submodels of class
#' \code{\link{summary.ctryvarest}}.
#'
#' @export
summary.bgvarest <- function(object, ci = .95, period = NULL, ctry = NULL, ...){
  
  if (is.null(ctry)) {
    pos <- 1:length(object)
  } else {
    pos <- which(names(object) %in% ctry) 
  }
  
  for (i in pos) {
    object[[i]] <- summary.ctryvarest(object[[i]], ci = ci, period = period, ...)
  }
  
  class(object) <- c("summary.bgvarest", "list")
  return(object)
}
