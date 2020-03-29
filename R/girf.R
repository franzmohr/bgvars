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
#' @examples 
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' global_data <- gvar2016$global_data
#' region_weights <- gvar2016$region_weights
#' weight_data <- gvar2016$weight_data
#' 
#' # Take first difference of all variables y and Dp
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"))
#' 
#' # Generate EA area region with 2 year, rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "MX", "GB"),
#'                                      p_domestic = 1,
#'                                      p_foreign = 1,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         model_specs = model_specs)
#' 
#' # Add priors
#' object <- add_priors(object)
#' 
#' # Estimate GVAR model
#' gvar_est <- estimate_gvar(object, iterations = 100, burnin = 10, thin = 2)
#' # Note that the number of iterations and burnin should be much higher.
#' 
#' # Solve GVAR
#' gvar_solved <- solve_gvar(gvar_est)
#' 
#' # GIRF
#' gvar_girf <- girf(gvar_solved, impulse = c("US", "r"), response = c("EA", "y"))
#' 
#' # Plot GIRF
#' plot(gvar_girf)
#' 
#' @export
girf <- function(object, impulse, response, n.ahead = 5,
                 shock = "sd", ci = .95, cumulative = FALSE, mc.cores = NULL) {
  
  if (!"bgvar" %in% class(object)) {
    stop("Object must be of class 'bgvar'.")
  }
  
  impulse_attr <- impulse
  response_attr <- response
  impulse <- which(object$index[,"country"] == impulse[1] & object$index[, "variable"] == impulse[2])
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(object$index[,"country"] == response[1] & object$index[, "variable"] == response[2])
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- ncol(object$data$X)
  store <- nrow(object$g0_i)
  
  a <- NULL
  for (i in 1:store) {
    a[[i]] <- list(g0_i = matrix(object$g0_i[i, ], k),
                   g = matrix(object$g[i, ], k),
                   sigma = matrix(object$sigma[i, ], k))
    if (shock == "sd") {
      a[[i]]$shock <- sqrt(a[[i]]$sigma[impulse, impulse])
    }
    if (shock == "nsd") {
      a[[i]]$shock <- -sqrt(a[[i]]$sigma[impulse, impulse])      
    }
  }
  
  if (is.null(mc.cores)) {
    result <- lapply(a, .ir, h = n.ahead, impulse = impulse, response = response, full = FALSE) 
  } else {
    result <- parallel::mclapply(a, .ir, h = n.ahead, impulse = impulse,
                                 response = response, full = FALSE, mc.cores = mc.cores)
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