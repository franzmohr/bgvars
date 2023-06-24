#' Impulse Response Function for a GVAR Submodel
#' 
#' Computes the impulse response coefficients of an object of class \code{"ctryvarest"} for
#' \code{n.ahead} steps.
#' 
#' @param x an object of class \code{"ctryvarest"}, usually, a result of a call to
#' \code{\link{draw_posterior.gvarsubmodels}} or \code{\link{ctryvec_to_ctryvar}}.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param shock size of the shock.
#' @param type type of the impulse response. Possible choices are forecast error \code{"feir"}
#' (default), orthogonalised \code{"oir"}, structural \code{"sir"}, generalised \code{"gir"},
#' and structural generalised \code{"sgir"} impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the posterior impulse response function. Defaults to \code{FALSE} so that
#' the median and the credible intervals of the posterior draws are returned.
#' @param period integer. Index of the period, for which the IR should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param ctry character. Name of the element in argument \code{x}, for which IRFs
#' should be obtained. If \code{NULL} (default), all submodels are used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A time-series object of class \code{"bvarirf"} and if \code{keep_draws = TRUE} a simple matrix.
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
irf.bgvarest <- function(x, impulse, response, n.ahead = 5, ci = .95, shock = 1,
                         type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL,
                         ctry = NULL, ...) {
  
  pos <- 1:length(x)
  if (!is.null(ctry)) {
    pos <- which(names(x) %in% ctry) 
  }
  
  if (length(pos) == 0) {
    stop("Specified countries not available.")
  }
  
  result <- list()
  for (i in pos) {
    
    # Skip tests if posterior simulation was not successful
    if (!is.null(x[[i]][["error"]])) {
      if (x[[i]][["error"]]) {
        next 
      }
    }
    
    result[[i]] <- irf(x[[i]], impulse = impulse, response = response, n.ahead = n.ahead, ci = ci,
                       shock = shock, type = type, cumulative = cumulative,
                       keep_draws = keep_draws, period = period, ...)
  }
  
  class(result) <- c("bgvarestirf", "list")
  return(result)
}