#' Posterior Simulation Country-Specific VARX Models of a GVAR Model
#' 
#' Produces draws from the posterior distributions of Bayesian VARX models.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' in combination with \code{\link{add_priors}}.
#' 
#' @details The function implements commonly used posterior simulation algorithms for Bayesian VAR models.
#' It can produce posterior draws for standard BVAR models with independent normal-Wishart priors, which can
#' be augmented by stochastic search variable selection (SSVS) as proposed by Geroge et al. (2008) or Bayesian
#' variable selection (BVS) as proposed in Korobilis (2013). Both SSVS and BVS can also be applied to the
#' covariances of the error term.
#' 
#' The implementation follows the description in Chan et al. (2019), George et al. (2008) and Korobilis (2013).
#' For all approaches the SUR form of a VAR model is used to obtain posterior draws. The algorithm is implemented
#' in C++ to reduce calculation time.
#' 
#' The function also supports structural BVAR models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VAR model, the structural coefficients are drawn jointly with the other coefficients.
#' 
#' @return An object of class \code{"ctryvarest"}.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \url{https://doi.org/10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \url{https://doi.org/10.1002/jae.1271}
#' 
#' @export
bvarxpost <- function(object) {
  
  if (object[["model"]][["tvp"]]) {
    object <- .bgvartvpalg(object)
  } else {
    object <- .bgvaralg(object)
  }
  
  k <- NCOL(object[["data"]][["Y"]])
  tt <- NROW(object[["data"]][["Y"]])
  n_dom <- 0
  if (!is.null(object$model$domestic)) {
    n_dom <- length(object$model$domestic[["variables"]]) * object$model$domestic[["lags"]]
  }
  n_for <- length(object$model$foreign[["variables"]]) * object$model$foreign[["lags"]]
  n_glo <- 0
  if (!is.null(object$model$global)) {
    n_glo <- length(object$model$global[["variables"]]) * object$model$global[["lags"]]
  }
  n_det <- 0
  if (!is.null(object$model$deterministic)) {
    n_det <- length(object$model$deterministic)
  }
  tot_pars <- n_dom + n_for + n_glo + n_det

  tvp_a0 <- FALSE
  tvp_dom <- FALSE
  tvp_for <- FALSE
  tvp_glo <- FALSE
  tvp_det <- FALSE
  tvp_sigma <- FALSE
  structural <- FALSE
  
  mc_spec <- c(object[["model"]][["burnin"]] + 1, object[["model"]][["iterations"]], 1)

  if (!is.null(object[["posteriors"]][["sigma"]][["lambda"]])) {
    sigma_lambda <- matrix(diag(NA_real_, k), k * k, object[["model"]][["iterations"]])
    sigma_lambda[which(lower.tri(diag(1, k))), ] <- object[["posteriors"]][["sigma"]][["lambda"]]
    sigma_lambda[which(upper.tri(diag(1, k))), ] <- object[["posteriors"]][["sigma"]][["lambda"]]
    object[["posteriors"]][["sigma"]][["lambda"]] <- sigma_lambda
    rm(sigma_lambda)
  }

  #### Combine coefficients ####
  
  A0 <- NULL
  if (object[["model"]][["structural"]]) {

    pos <- which(lower.tri(diag(1, k)))
    draws <- object[["model"]][["iterations"]]

    if (is.list(object[["posteriors"]][["a0"]])) {

      if ("coeffs" %in% names(object[["posteriors"]][["a0"]])) {
        if (object[["model"]][["tvp"]]) {
          A0[["coeffs"]] <- matrix(diag(1, k), k * k * tt, draws)
          A0[["coeffs"]][rep(0:(tt - 1) * k * k, each = length(pos)) + rep(pos, tt), ] <- object[["posteriors"]][["a0"]][["coeffs"]]
        } else {
          A0[["coeffs"]] <- matrix(diag(1, k), k * k, draws)
          A0[["coeffs"]][rep(0:(draws - 1) * k * k, each = length(pos)) + pos ] <- object[["posteriors"]][["a0"]][["coeffs"]]
        }
      }

      if ("sigma" %in% names(object[["posteriors"]][["a0"]])) {
        A0[["sigma"]] <- matrix(0, k * k, draws)
        A0[["sigma"]][pos, ] <- object[["posteriors"]][["a0"]][["sigma"]]
      }

      if ("lambda" %in% names(object[["posteriors"]][["a0"]])) {
        A0[["lambda"]] <- matrix(diag(1, k), k * k, draws)
        A0[["lambda"]][pos, ] <- object[["posteriors"]][["a0"]][["lambda"]]
        A0[["lambda"]][-pos, ] <- NA_real_
      }

    } else {
      A0 <- matrix(diag(1, k), k * k, draws)
      A0[pos, ] <- object[["posteriors"]][["a0"]]
    }

    object[["posteriors"]][["a0"]] <- NULL

    if(!is.null(A0)) {
      if (is.list(A0)) {
        if ("coeffs" %in% names(A0)) {
          n_a0 <- nrow(A0[["coeffs"]])
        }
      } else {
        n_a0 <- nrow(A0)
      }
      if (n_a0 / tt >= 1) {
        tvp_a0 <- TRUE
        n_a0 <- n_a0 / tt
      }
      if (n_a0 %% (k * k) != 0) {
        stop("Row number of coefficient draws of 'A0' is not k^2 or multiples thereof.")
      }
      structural <- TRUE

      object$posteriors <- c(object$posteriors, .country_fill_helper(A0, tvp_a0, n_a0, tt, "a0"))
    }
  }

  if(!is.null(object$posteriors$domestic)) {
    if (is.list(object$posteriors$domestic)) {
      if ("coeffs" %in% names(object$posteriors$domestic)) {
        n_dom <- nrow(object$posteriors$domestic[["coeffs"]])
      }
    } else {
      n_dom <- nrow(object$posteriors$domestic)
    }
    if ((n_dom / tt) %% k^2 == 0 & n_dom / tt >= 1) {
      tvp_dom <- TRUE
      n_dom <- n_dom / tt
    }

    domestic <- object$posteriors$domestic
    object$posteriors$domestic <- NULL
    object$posteriors <- c(object$posteriors, .country_fill_helper(domestic, tvp_dom, n_dom, tt, "domestic"))
    rm(domestic)
  }

  if(!is.null(object$posteriors$foreign)) {
   if (is.list(object$posteriors$foreign)) {
      if ("coeffs" %in% names(object$posteriors$foreign)) {
        n_for <- nrow(object$posteriors$foreign[["coeffs"]])
      }
    } else {
      n_for <- nrow(object$posteriors$foreign)
    }
    if ((n_for / tt) %% k == 0) {
      tvp_for <- TRUE
      n_for <- n_for / tt
    }

    foreign <- object$posteriors$foreign
    object$posteriors$foreign <- NULL
    object$posteriors <- c(object$posteriors, .country_fill_helper(foreign, tvp_for, n_for, tt, "foreign"))
    rm(foreign)
  }

  if(!is.null(object$posteriors$global)) {
    if (is.list(object$posteriors$global)) {
      if ("coeffs" %in% names(object$posteriors$global)) {
        n_glo <- nrow(object$posteriors$global[["coeffs"]])
      }
    } else {
      n_glo <- nrow(object$posteriors$global)
    }
    if ((n_glo / tt) %% k == 0) {
      tvp_glo <- TRUE
      n_glo <- n_glo / tt
    }

    global <- object$posteriors$global
    object$posteriors$global <- NULL
    object$posteriors <- c(object$posteriors, .country_fill_helper(global, tvp_glo, n_glo, tt, "global"))
    rm(global)
  }

  if(!is.null(object$posteriors$deterministic)) {

    if (is.list(object$posteriors$deterministic)) {
      if ("coeffs" %in% names(object$posteriors$deterministic)) {
        n_det <- NROW(object$posteriors$deterministic[["coeffs"]])
      }
    } else {
      n_det <- NROW(object$posteriors$deterministic)
    }
    if ((n_det / tt) %% k == 0 & n_det / tt >= 1) {
      tvp_det <- TRUE
      n_det <- n_det / tt
    }

    deterministic <- object$posteriors$deterministic
    object$posteriors$deterministic <- NULL
    object$posteriors <- c(object$posteriors, .country_fill_helper(deterministic, tvp_det, n_det, tt, "deterministic"))
    rm(deterministic)
  }

  if(!is.null(object$posteriors$sigma)) {
    if (is.list(object$posteriors$sigma)) {
      if ("coeffs" %in% names(object$posteriors$sigma)) {
        n_sigma <- nrow(object$posteriors$sigma[["coeffs"]])
      }
    } else {
      n_sigma <- nrow(object$posteriors$sigma)
    }
    if ((n_sigma / tt) %% k == 0 & n_sigma / tt >= 1) {
      tvp_sigma <- TRUE
      n_sigma <- n_sigma / tt
    }
    if (n_sigma %% (k * k) != 0) {
      stop("Row number of coefficient draws of 'Sigma' is not k^2 or multiples thereof.")
    }

    sigma <- object$posteriors$sigma
    object$posteriors$sigma <- NULL
    object$posteriors <- c(object$posteriors, .country_fill_helper(sigma, tvp_sigma, n_sigma, tt, "sigma"))
    rm(sigma)
  }
  
  vars <- c("a0",
            "domestic",
            "foreign",
            "global",
            "determinisitc",
            "sigma")
  
  for (i in vars) {
    if (!is.null(object[["posteriors"]][[i]])) {
      attr(object[["posteriors"]][[i]], "mcpar") <- mc_spec
    }
  }

  class(object) <- append("ctryvarest", class(object))
  
  return(object)
}