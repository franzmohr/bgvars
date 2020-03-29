#' Generalised Forecast Error Variance Decomposition
#' 
#' Produces the generalised forecast error variance decomposition of a Bayesian GVAR model.
#' 
#' @param object an object of class \code{"bgvar"}, usually, a result of a call to \code{\link{solve_gvar}}.
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
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
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
#' # Generate EA area region with 3 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 3,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 3 year rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 3,
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
#' # GFEVD
#' gvar_gfevd <- gfevd(gvar_solved, response = c("EA", "y"))
#'
#' # Plot GFEVD
#' plot(gvar_gfevd)
#' 
#' @export
gfevd <- function(object, response, n.ahead = 5, normalise_gir = FALSE, mc.cores = NULL) {
  if (!"bgvar" %in% class(object)) {
    stop("Object must be of class 'bgvar'.")
  }
  if (is.null(object$sigma)) {
    stop("The 'bgvar' object must include draws of the variance-covariance matrix Sigma.")
  }
  
  response <- which(object$index[, "country"] == response[1] & object$index[, "variable"] == response[2])
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- sqrt(NCOL(object$g0_i))
  store <- NROW(object$g0_i)
  
  # Produce FEIR
  a <- NULL # Prepare data for lapply
  for (i in 1:store) {
    a[[i]] <- list(g0_i = matrix(object$g0_i[i, ], k),
                   g = matrix(object$g[i, ], k),
                   sigma = matrix(object$sigma[i, ], k))
    a[[i]]$shock <- 0
  }
  if (is.null(mc.cores)) {
    phi <- lapply(a, .ir, h = n.ahead, impulse = 0, response = 0, full = TRUE)
  } else {
    phi <- parallel::mclapply(a, .ir, h = n.ahead, impulse = 0, response = 0,
                              full = FALSE, mc.cores = mc.cores)
  }

  # Produce GFEVD
  ej_t <- matrix(0, 1, k)
  ej_t[, response] <- 1
  result <- matrix(NA, (n.ahead + 1) * k, store)
  for (j in 1:store) {
    P <- a[[j]]$g0_i %*% a[[j]]$sigma
    numerator <- matrix(NA, n.ahead + 1, k)
    mse <- matrix(NA, n.ahead + 1, 1)
    numerator[1, ] <- (ej_t %*% phi[[j]][1:k, ] %*% P)^2
    mse[1,] <- ej_t %*% phi[[j]][1:k,] %*% tcrossprod(P, a[[j]]$g0_i) %*% t(phi[[j]][1:k,]) %*% t(ej_t)
    for (i in 2:(n.ahead + 1)) {
      numerator[i, ] <- numerator[i - 1, ] + (ej_t %*% phi[[j]][(i - 1) * k + 1:k, ] %*% P)^2
      mse[i, ] <- mse[i - 1,] + ej_t %*% phi[[j]][(i - 1) * k + 1:k,] %*% tcrossprod(P, a[[j]]$g0_i) %*% t(phi[[j]][(i - 1) * k + 1:k,]) %*% t(ej_t)
    }
    numerator <- numerator / a[[j]]$sigma[response, response]
    result[, j] <- matrix(numerator / matrix(mse, n.ahead + 1, k))
  }
  # Get means
  result <- matrix(apply(result, 1, mean), n.ahead + 1)
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