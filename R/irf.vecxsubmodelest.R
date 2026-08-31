#' Impulse Response Function for a GVAR Submodel
#' 
#' Computes the impulse response coefficients of an object of class 'vecxsubmodelest' for
#' \code{n_ahead} steps.
#' 
#' @param x an object of class 'vecxsubmodelest', usually, a result of a call to
#' \code{\link{draw_posterior.gvarsubmodels}} or \code{\link[bvartools]{bvec_to_bvar}}.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n_ahead number of steps ahead.
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
irf.vecxsubmodelest <- function(x, impulse, response, n_ahead = 5, ci = .95, shock = 1,
                                type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL, ...) {
  
  # Transform to bvar-object ----
  k_domestic <- x[["model"]][["k_domestic"]]
  p_domestic <- x[["model"]][["p_domestic"]]
  k_foreign <- x[["model"]][["k_foreign"]]
  p_foreign <- x[["model"]][["p_foreign"]]
  m <- x[["model"]][["m"]]
  s <- x[["model"]][["s"]]
  n_gamma_domestic <- k_domestic * k_domestic * (p_domestic - 1)
  r <- x[["model"]][["rank"]]
  n_alpha <- r * k_domestic
  tt <- nrow(x[["data"]][["y"]])
  draws <- x[["model"]][["iterations"]]
  tvp <- x[["model"]][["tvp"]]
  n_z <- x[["data"]][["z"]]
  
  # Modelled variables ----
  draws_gamma <- NULL
  if (p_domestic > 1) {
    pos_a <- n_alpha + 1:n_gamma_domestic
    if (tvp) {
      pos_a <- rep(pos_a, tt) + rep(0:(tt - 1), each = length(pos_a)) * n_z
    }
    draws_gamma[["coeffs"]] <- t(x[["posteriors"]][["a"]][["coeffs"]][, pos_a])
  }
  
  # Structural ----
  draws_a0 <- NULL
  if (x[["model"]][["structural"]]) {
    stop("Structural models not implemented yet.")
    pos <- which(lower.tri(diag(1, k_domestic)))
    pos_a0 <- (n_z - length(pos) + 1):n_z
    
    if (x[["model"]][["tvp"]]) {
      stop("implement tvp")
      pos_a0_long <- rep(pos_a0, tt) + rep(0:(tt - 1), each = length(pos_a0)) * n_z
      draws_a0[["coeffs"]] <- matrix(diag(1, k), k * k * tt, draws)
      draws_a0[["coeffs"]][rep(0:(tt - 1) * k * k, each = length(pos)) + rep(pos, tt), ] <- t(x[["posteriors"]][["a"]][["coeffs"]][, pos_a0_long])
    } else {
      draws_a0[["coeffs"]] <- matrix(diag(1, k_domestic), k_domestic * k_domestic, draws)
      draws_a0[["coeffs"]][pos, ] <- t(x[["posteriors"]][["a"]][["coeffs"]][, pos_a0])
    }
  }
  
  
  # Sigma----
  draws_sigma <- NULL
  draws_sigma[["coeffs"]] <- t(x[["posteriors"]][["sigma"]][["coeffs"]])
  
  # Create bvar object
  bvar_object <- bvar(data = NULL,
                      exogen = NULL,
                      y = x[["data"]][["y"]],
                      x = x[["data"]][["x"]],
                      z = x[["data"]][["z"]],
                      A0 = draws_a0,
                      A = draws_gamma,
                      Sigma = draws_sigma)
  
  result <- irf(bvar_object, impulse = impulse, response = response,
                n_ahead = n_ahead, ci = ci, shock = shock, type = type,
                cumulative = cumulative, keep_draws = keep_draws, period = period, ...)
  
  
  return(result)
}