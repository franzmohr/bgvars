#' Time Series Windows
#' 
#' Forwards model input to the same function for individual models.
#' 
#' @param x an object of class 'gvarmodel'.
#' @param start the start time of the period of interest.
#' @param end the end time of the period of interest.
#' 
#' @return An object of class 'gvarmodel'.
#' 
#' @export
#' @method window gvarmodel
window.gvarmodel <- function(x, start = NULL, end = NULL, ...) {
  
  for (i in 1:length(x[["submodels"]])) {
    x[["submodels"]][[i]] <- stats::window(x = x[["submodels"]][[i]], start = start, end = end, ...)
  }
  
  return(x)
}