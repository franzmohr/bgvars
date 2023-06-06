#' Summarising Country-Specific VECX Submodels of a GVAR Model
#'
#' summary method for class \code{"bgvecest"}.
#'
#' @param object an object of class \code{"bgvecest"}, usually, a result of a call to
#' \code{\link{draw_posterior}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param period integer. Index of the period, for which the summary statistics should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param ctry character. Name of the element in argument \code{object}, for which summary
#' statistics should be obtained. If \code{NULL} (default), all country submodels are used.
#' @param x an object of class \code{"summary.bgvecest"}, usually, a result of a call to
#' \code{\link{summary.bgvecest}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.bgvecest} returns a list of class \code{"summary.bgvecest"},
#' which contains summary statistics of country-specific submodels of class
#' \code{\link{summary.ctryvecest}}.
#'
#' @export
summary.bgvecest <- function(object, ci = .95, period = NULL, ctry = NULL, ...){
  
  pos <- 1:length(object)
  if (!is.null(ctry)) {
    pos <- which(names(object) %in% ctry) 
  }
  
  if (length(pos) == 0) {
    stop("Specified countries not available.")
  }
  
  result <- NULL
  for (i in pos) {
    temp <- summary.ctryvecest(object[[i]], ci = ci, period = period, ...) 
    result <- c(result, list(temp))
    rm(temp)
  }
  names(result) <- names(object)[pos]
  
  class(result) <- c("summary.bgvecest", "list")
  return(result)
}
