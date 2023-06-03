#' Gernalised Impulse Response Function
#' 
#' Computes the generalised impulse response coefficients of an object of class \code{"bgvar"} for
#' `n.ahead` steps.
#' 
#' @param object an object of class \code{"bgvar"}, usually, a result of a call to \code{\link{solve_gvar}}.
#' @param impulse a character vector of the impulse country and variable, respectively.
#' @param response a character vector of the response country and variable, respectively.y.
#' @param n.ahead number of steps ahead.
#' @param shock size of the shock. Either \code{"sd"} (default) for a standard deviation of
#' the impulse variable, \code{"nsd"} for a negativ standard deviation, or a numeric.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' 
#' @return A time-series object of class \code{"bgvarirf"}.
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' 
#' 
#' @export
girf <- function(object, impulse, response, n.ahead = 5,
                 shock = "sd", ci = .95, cumulative = FALSE, mc.cores = NULL) {
  
  # rm(list = ls()[-which(ls() == "object")]); impulse = c("US", "r"); response = c("US", "y"); n.ahead = 10; shock = "sd"; ci = .95; cumulative = FALSE; mc.cores = NULL
  
  if (!"bgvar" %in% class(object)) {
    stop("Object must be of class 'bgvar'.")
  }
  
  # Identify position of impulse and response in global matrix
  impulse_attr <- impulse
  response_attr <- response
  impulse <- which(object$index[,"country"] == impulse[1] & object$index[, "variable"] == impulse[2])
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(object$index[,"country"] == response[1] & object$index[, "variable"] == response[2])
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- ncol(object$data$endogen)
  store <- nrow(object$a0)
  
  # Generate a data object that can be passed to .ir
  a <- NULL
  for (i in 1:store) {
    a[[i]] <- list(a0 = matrix(object$a0[i, ], k),
                   a = matrix(object$a[i, ], k),
                   sigma = matrix(object$sigma[i, ], k))
    if (class(shock) == "character") {
      if (shock == "sd") {
        a[[i]]$shock <- sqrt(a[[i]]$sigma[impulse, impulse])
      }
      if (shock == "nsd") {
        a[[i]]$shock <- -sqrt(a[[i]]$sigma[impulse, impulse])      
      } 
    } else {
      a[[i]]$shock <- shock
    }
  }
  
  # Generate impulse responses
  if (is.null(mc.cores)) {
    result <- lapply(a, .ir, h = n.ahead, impulse = impulse, response = response) 
  } else {
    result <- parallel::mclapply(a, .ir, h = n.ahead, impulse = impulse,
                                 response = response, mc.cores = mc.cores)
  }
  
  result <- t(matrix(unlist(result), n.ahead + 1))
  
  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  pr <- c(ci_low, .5, ci_high)
  result <- stats::ts(t(apply(result, 2, stats::quantile, probs = pr)))
  stats::tsp(result) <- c(0, n.ahead, stats::tsp(result)[3])
  
  attr(result, "impulse") <- impulse_attr
  attr(result, "response") <- response_attr
  
  class(result) <- append("bgvarirf", class(result))
  return(result)
}