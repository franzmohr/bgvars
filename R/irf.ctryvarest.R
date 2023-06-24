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
irf.ctryvarest <- function(x, impulse, response, n.ahead = 5, ci = .95, shock = 1,
                           type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL, ...) {

  if (!type %in% c("feir", "oir", "gir", "sir", "sgir")) {
    stop("Argument 'type' not known.")
  }
  
  if (!"ctryvarest" %in% class(x)) {
    stop("Argument 'x' must be of class 'ctryvarest'.")
  }
  
  if (is.null(x[["posteriors"]][["domestic"]]) & !type %in% c("sir", "sgir")) {
    stop("Impulse responses only supported for models with p > 0, i.e. argument 'x' must contain element 'x$posteriors$domestic' or structural models.")
  }
  
  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (is.null(x[["posteriors"]][["a0"]])) {
      stop("Structural IR requires that draws of 'a0' are available.")
    }
    need_A0 <- TRUE
  }
  
  if (!(is.numeric(shock) | shock %in% c("sd", "nsd"))) {
    stop("Invalid specification of argument 'shock'.")
  }
  
  if (type %in% c("oir", "gir", "sgir") | shock %in% c("sd", "nsd")) {
    if (is.null(x[["posteriors"]][["sigma"]])) {
      stop("OIR, GIR, SGIR require that the 'bvar' x contains draws of 'Sigma'.")
    }
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
  }
  
  impulse_attr <- impulse
  response_attr <- response
  impulse <- which(dimnames(x[["data"]][["domestic"]])[[2]] == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(dimnames(x[["data"]][["domestic"]])[[2]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k_dom <- NCOL(x[["data"]][["Y"]])
  tt <- NROW(x[["data"]][["Y"]])
  tvp <- x[["model"]][["tvp"]]
  if (any(unlist(tvp))) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
  }
  
  store <- NA
  vars <- c("a0", "domestic", "foreign", "global", "deterministic", "sigma")
  for (i in vars) {
    if (is.na(store)) {
      if (!is.null(x[["posteriors"]][[i]])) {
        if (tvp) {
          store <- nrow(x[["posteriors"]][[i]][[1]])
        } else {
          store <- nrow(x[["posteriors"]][[i]]) 
        }
      }   
    }
  }
  
  A <- NULL
  for (i in 1:store) {
    temp <- NULL
    if (!is.null(x[["posteriors"]][["domestic"]])) {
      if (tvp) {
        temp[["A"]] <- matrix(x[["posteriors"]][["domestic"]][[period]][i, ], k_dom)
      } else {
        temp[["A"]] <- matrix(x[["posteriors"]][["domestic"]][i, ], k_dom) 
      }
    } else {
      temp[["A"]] <- matrix(0, k, k)
    }
    if (need_Sigma) {
      if (x[["model"]][["sv"]]) {
        temp[["Sigma"]] <- matrix(x[["posteriors"]][["sigma"]][[period]][i, ], k_dom)
      } else {
        temp[["Sigma"]] <- matrix(x[["posteriors"]][["sigma"]][i, ], k_dom) 
      }
    }
    
    # Shock
    if (is.numeric(shock)) {
        temp[["shock"]] <- shock
    } else {
      if (type == "oir") {
        temp[["shock"]] <- diag(chol(temp[["posteriors"]][["sigma"]]))[impulse]
      } else {
        temp[["shock"]] <- sqrt(diag(temp[["posteriors"]][["sigma"]])[impulse]) 
      }
      
      if (shock == "nsd") {
        temp[["shock"]] <- -temp[["shock"]]
      } 
    }
      
    if (need_A0) {
      if (tvp) {
        temp[["A0"]] <- matrix(x[["posteriors"]][["a0"]][[period]][i, ], k_dom)
      } else {
        temp[["A0"]] <- matrix(x[["posteriors"]][["a0"]][i, ], k_dom) 
      }
    }
    
    A[[i]] <- temp
  }

  result <- lapply(A, .ir, h = n.ahead, type = type,
                   impulse = impulse, response = response)
  
  result <- t(matrix(unlist(result), n.ahead + 1))
  
  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }
  
  if (!keep_draws) {
    ci_low <- (1 - ci) / 2
    ci_high <- 1 - ci_low
    pr <- c(ci_low, .5, ci_high)
    result <- stats::ts(t(apply(result, 2, stats::quantile, probs = pr)), start = 0, frequency = 1) 
  }
  
  attr(result, "impulse") <- impulse_attr
  attr(result, "response") <- response_attr
  if (!is.null(x[["model"]][["rank"]])) {
    attr(result, "rank") <- x[["model"]][["rank"]]    
  }
  
  class(result) <- append("bvarirf", class(result))
  return(result)
}