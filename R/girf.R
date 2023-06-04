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
#' # Load data
#' data("gvar2016")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2016$country_data,
#'                        weight_data = gvar2016$weight_data,
#'                        region_weights = gvar2016$region_weights,
#'                        regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")), period = 3)
#'
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2016$global_data
#' 
#' # Difference series to make them stationary
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
#' global_data <- diff_variables(global_data, multi = 100)
#' 
#' # Create weights
#' weight_data <- create_weights(weight_data, period = 3, country_data = country_data)
#' 
#' # Create general model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      domestic = list(variables = c("y", "Dp", "r"), lags = 1:2),
#'                                      foreign = list(variables = c("y", "Dp", "r"), lags = 1:2),
#'                                      global = list(variables = "poil", lags = 1:2),
#'                                      deterministic = list("const" = TRUE),
#'                                      countries = c("EA", "US", "GB", "CA", "JP", "IN"),
#'                                      type = "VAR",
#'                                      tvp = FALSE, sv = FALSE,
#'                                      iterations = 10,
#'                                      burnin = 10)
#' # Note that the number of iterations and burnin should be much higher!
#'                                      
#' # Country-specific specifications
#' model_specs[["US"]][["domestic"]][["variables"]] <- c("y", "Dp", "r")
#' model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")
#' 
#' # Create all country models
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#'
#' 
#' # Add priors
#' temp_model <- add_priors(country_models,
#'                          coef = list(v_i = 1 / 9, v_i_det = 1 / 10),
#'                          sigma = list(df = 3, scale = 1))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(temp_model)
#' 
#' # Solve GVAR
#' gvar <- combine_submodels(object)
#' 
#' # Obtain GIRF
#' gvar_girf <- girf(gvar, impulse = c("EA", "r"), response = c("EA", "y"), n.ahead = 20)
#' 
#' # Plot GIRF
#' plot(gvar_girf)
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