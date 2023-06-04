#' Generalised Forecast Error Variance Decomposition
#' 
#' Produces the generalised forecast error variance decomposition of a Bayesian GVAR model.
#' 
#' @param object an object of class \code{"bgvar"}, usually, a result of a call to \code{\link{combine_submodels}}.
#' @param response a character vector of the response country and variable, respectively.
#' @param n.ahead number of steps ahead.
#' @param normalise_gir logical. Should the GFEVD be normalised?
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' 
#' @details For the global VAR model
#' \deqn{y_t = \sum_{l = 1}^{p} G_l y_{t - j} + G^{-1}_{0} u_t}
#' with \eqn{u_t \sim \Sigma} and \eqn{G_i} as \eqn{K \times K} coefficient matrices
#' the function produces the generalised structural forecast error variance decomposition as
#' \deqn{\omega^{GIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i G_0^{-1} \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i G_0^{-1} \Sigma G_0^{-1 \prime} \Phi_i^{\prime} e_j )},}
#' where \eqn{\Phi_i} is the forecast error impulse response for the \eqn{i}th period, \eqn{\Sigma}
#' is the variance-covariance matrix of the error term,
#' \eqn{e_j} is a selection vector for the response variable,
#' \eqn{e_k} is a selection vector for the impulse variable,
#' and \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#' 
#' Since GIR-based FEVDs do not add up to unity, they can be normalised by setting \code{normalise_gir = TRUE}.
#' 
#' @return A time-series object of class "bgvarfevd".
#' 
#' @examples
#' # Load data
#' data("gvar2016")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2016$country_data,
#'              weight_data = gvar2016$weight_data,
#'              region_weights = gvar2016$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2016$global_data
#' 
#' # Difference series to make them stationary
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
#' global_data <- diff_variables(global_data, multi = 100)
#' 
#' # Create time varying weights
#' weight_data <- create_weights(weight_data, period = 3, country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(
#'                  country_data = country_data,
#'                  global_data = global_data,
#'                  countries = c("US", "JP", "CA", "NO", "GB", "EA"), 
#'                  domestic = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  global = list(variables = c("poil"), lags = 1),
#'                  deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                  iterations = 10,
#'                  burnin = 10)
#' # Note that the number of iterations and burnin draws should be much higher!
#'                                      
#' # Overwrite country-specific specifications
#' model_specs[["US"]][["domestic"]][["variables"]] <- c("y", "Dp", "r")
#' model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")
#' 
#' # Create estimation objects
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#' 
#' # Add priors
#' models_with_priors <- add_priors(country_models,
#'                                  coef = list(v_i = 1 / 9, v_i_det = 1 / 10),
#'                                  sigma = list(df = 3, scale = .0001))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(models_with_priors)
#' 
#' # Solve GVAR
#' gvar <- combine_submodels(object)
#' 
#' # Obtain forecasts
#' gvar_fevd <- gfevd(gvar, response = c("US", "y"))
#' 
#' # Plot forecast
#' plot(gvar_fevd) 
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#'  
#' @export
gfevd <- function(object, response, n.ahead = 5, normalise_gir = FALSE, mc.cores = NULL) {
  
  # rm(list = ls()[-which(ls() == "object")]); response = c("US", "y"); n.ahead = 5; normalise_gir = FALSE; mc.cores = NULL
  
  if (!"bgvar" %in% class(object)) {
    stop("Object must be of class 'bgvar'.")
  }
  if (is.null(object$sigma)) {
    stop("The 'bgvar' object must include draws of the variance-covariance matrix Sigma.")
  }
  
  response <- which(object$index[, "country"] == response[1] & object$index[, "variable"] == response[2])
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- sqrt(NCOL(object$a0)) # Number of endogenous variables
  store <- NROW(object$a0) # Number of draws
  
  # Produce FEIR
  a <- NULL # Prepare data for lapply
  for (i in 1:store) {
    a[[i]] <- list(a0 = matrix(object$a0[i, ], k),
                   a = matrix(object$a[i, ], k),
                   sigma = matrix(object$sigma[i, ], k))
  }
  
  if (is.null(mc.cores)) {
    phi <- lapply(a, .vardecomp, h = n.ahead, response = response)
  } else {
    phi <- parallel::mclapply(a, .vardecomp, h = n.ahead, response = response,
                              mc.cores = mc.cores)
  }
  
  result <- matrix(rowMeans(matrix(unlist(phi), (n.ahead + 1) * k)), n.ahead + 1)

  # Normalise
  if (normalise_gir) {
    result <- t(apply(result, 1, function(x) {x / sum(x)}))
  }
  # Name columns
  dimnames(result) <- list(NULL, paste(object$index[, "country"], object$index[, "variable"], sep = "_"))
  # Turn into time-series object
  result <- stats::ts(result, start = 0, frequency = 1)
  # Define class
  class(result) <- append("bgvarfevd", class(result))
  return(result)
}