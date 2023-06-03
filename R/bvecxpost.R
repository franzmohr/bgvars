#' Posterior Simulation Country-Specific VECX Models of a GVAR Model
#' 
#' Produces draws from the posterior distributions of Bayesian VECX models.
#' 
#' @param object an object of class \code{"gvecsubmodels"}, usually, a result of a call to \code{\link{create_models}}
#' in combination with \code{\link{add_priors}}.
#' 
#' @details The function implements posterior simulation algorithms proposed in Koop et al. (2010)
#' and Koop et al. (2011), which place identifying restrictions on the cointegration space.
#' Both algorithms are able to employ Bayesian variable selection (BVS) as proposed in Korobilis (2013).
#' The algorithm of Koop et al. (2010) is also able to employ stochastic search variable selection (SSVS)
#' as proposed by Geroge et al. (2008).
#' Both SSVS and BVS can also be applied to the covariances of the error term. However, the algorithms
#' cannot be applied to cointegration related coefficients, i.e. to the loading matrix \eqn{\alpha} or
#' the cointegration matrix \eqn{beta}.
#' 
#' The implementation primarily follows the description in Koop et al. (2010). Chan et al. (2019),
#' George et al. (2008) and Korobilis (2013) were used to implement the variable selection algorithms.
#' For all approaches the SUR form of a VEC model is used to obtain posterior draws. The algorithm is implemented
#' in C++ to reduce calculation time.
#' 
#' The function also supports structural BVEC models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VEC model, the structural coefficients are drawn jointly with the other coefficients. No identifying
#' restrictions are made regarding the cointegration matrix.
#' 
#' @return An object of class \code{"bvar"}.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @export
bvecxpost <- function(object) {
  
  if (object[["model"]][["tvp"]]) {
    #object <- .bvecxtvpalg(object) # Use C++ code to draw posteriors 
  } else {
    object <- .bgvecalg(object)
  }
  
  k <- NCOL(object[["data"]][["Y"]])
  tt <- NROW(object[["data"]][["Y"]])
  
  tvp_a0 <- FALSE
  tvp_alpha <- FALSE
  tvp_beta_dom <- FALSE
  tvp_beta_for <- FALSE
  tvp_beta_glo <- FALSE
  tvp_beta_det <- FALSE
  tvp_gamma_dom <- FALSE
  tvp_gamma_for <- FALSE
  tvp_gamma_glo <- FALSE
  tvp_gamma_det <- FALSE
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
      
      object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(A0, tvp_a0, n_a0, tt, "a0"))
    }
  }
  
  # alpha ----
  if(!is.null(object[["posteriors"]][["alpha"]])) {
    
    if ("coeffs" %in% names(object[["posteriors"]][["alpha"]])) {
      n_alpha <- nrow(object[["posteriors"]][["alpha"]][["coeffs"]])
    }
    if ((n_alpha / tt) %% k^2 == 0 & n_alpha / tt >= 1) {
      tvp_alpha <- TRUE
      n_alpha <- n_alpha / tt
    }
    
    alpha <- object[["posteriors"]][["alpha"]]
    object[["posteriors"]][["alpha"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(alpha, tvp_alpha, n_alpha, tt, "alpha"))
    rm(alpha)
  }
  
  # Beta domestic ----
  if(!is.null(object[["posteriors"]][["beta_dom"]])) {
    
    if ("coeffs" %in% names(object[["posteriors"]][["beta_dom"]])) {
      n_beta_dom <- nrow(object[["posteriors"]][["beta_dom"]][["coeffs"]])
    }
    if ((n_beta_dom / tt) %% k^2 == 0 & n_beta_dom / tt >= 1) {
      tvp_beta_dom <- TRUE
      n_beta_dom <- n_beta_dom / tt
    }
    
    beta_dom <- object[["posteriors"]][["beta_dom"]]
    object[["posteriors"]][["beta_dom"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(beta_dom, tvp_beta_dom, n_beta_dom, tt, "beta_domestic"))
    rm(beta_dom)
  } else {
    object[["posteriors"]][["beta_dom"]] <- NULL
  }
  
  # Beta foreign ----
  if(!is.null(object[["posteriors"]][["beta_for"]])) {
    
    if ("coeffs" %in% names(object[["posteriors"]][["beta_for"]])) {
      n_beta_for <- nrow(object[["posteriors"]][["beta_for"]][["coeffs"]])
    }
    if ((n_beta_for / tt) %% k^2 == 0 & n_beta_for / tt >= 1) {
      tvp_beta_for <- TRUE
      n_beta_for <- n_beta_for / tt
    }
    
    beta_for <- object[["posteriors"]][["beta_for"]]
    object[["posteriors"]][["beta_for"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(beta_for, tvp_beta_for, n_beta_for, tt, "beta_foreign"))
    rm(beta_for)
  } else {
    object[["posteriors"]][["beta_for"]] <- NULL
  }
  
  # Beta global ----
  if(!is.null(object[["posteriors"]][["beta_glo"]])) {
    
    if ("coeffs" %in% names(object[["posteriors"]][["beta_glo"]])) {
      n_beta_glo <- nrow(object[["posteriors"]][["beta_glo"]][["coeffs"]])
    }
    if ((n_beta_glo / tt) %% k^2 == 0 & n_beta_glo / tt >= 1) {
      tvp_beta_glo <- TRUE
      n_beta_glo <- n_beta_glo / tt
    }
    
    beta_glo <- object[["posteriors"]][["beta_glo"]]
    object[["posteriors"]][["beta_glo"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(beta_glo, tvp_beta_glo, n_beta_glo, tt, "beta_global"))
    rm(beta_glo)
  } else {
    object[["posteriors"]][["beta_glo"]] <- NULL
  }
  
  # Beta deterministic ----
  if(!is.null(object[["posteriors"]][["beta_d"]])) {
    
    if ("coeffs" %in% names(object[["posteriors"]][["beta_d"]])) {
      n_beta_det <- nrow(object[["posteriors"]][["beta_d"]][["coeffs"]])
    }
    if ((n_beta_det / tt) %% k^2 == 0 & n_beta_det / tt >= 1) {
      tvp_beta_det <- TRUE
      n_beta_det <- n_beta_det / tt
    }
    
    beta_d <- object[["posteriors"]][["beta_d"]]
    object[["posteriors"]][["beta_d"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(beta_d, tvp_beta_det, n_beta_det, tt, "beta_deterministic"))
    rm(beta_d)
  } else {
    object[["posteriors"]][["beta_d"]] <- NULL
  }
  
  # Gamma domestic ----
  if(!is.null(object[["posteriors"]][["gamma_dom"]])) {
    if (is.list(object[["posteriors"]][["gamma_dom"]])) {
      if ("coeffs" %in% names(object[["posteriors"]][["gamma_dom"]])) {
        n_gamma_dom <- nrow(object[["posteriors"]][["gamma_dom"]][["coeffs"]])
      }
    } else {
      n_gamma_dom <- nrow(object[["posteriors"]][["gamma_dom"]])
    }
    if ((n_gamma_dom / tt) %% k^2 == 0 & n_gamma_dom / tt >= 1) {
      tvp_gamma_dom <- TRUE
      n_gamma_dom <- n_gamma_dom / tt
    }
    
    gamma_dom <- object[["posteriors"]][["gamma_dom"]]
    object[["posteriors"]][["gamma_dom"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(gamma_dom, tvp_gamma_dom, n_gamma_dom, tt, "gamma_domestic"))
    rm(gamma_dom)
  } else {
    object[["posteriors"]][["gamma_dom"]] <- NULL
  }
  
  # Gamma foreign ----
  if(!is.null(object[["posteriors"]][["gamma_for"]])) {
    if (is.list(object[["posteriors"]][["gamma_for"]])) {
      if ("coeffs" %in% names(object[["posteriors"]][["gamma_for"]])) {
        n_gamma_for <- nrow(object[["posteriors"]][["gamma_for"]][["coeffs"]])
      }
    } else {
      n_gamma_for <- nrow(object[["posteriors"]][["gamma_for"]])
    }
    if ((n_gamma_for / tt) %% k^2 == 0 & n_gamma_for / tt >= 1) {
      tvp_gamma_for <- TRUE
      n_gamma_for <- n_gamma_for / tt
    }
    
    gamma_for <- object[["posteriors"]][["gamma_for"]]
    object[["posteriors"]][["gamma_for"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(gamma_for, tvp_gamma_for, n_gamma_for, tt, "gamma_foreign"))
    rm(gamma_for)
  }
  
  # Gamma global ----
  if(!is.null(object[["posteriors"]][["upsilon"]])) {
    if (is.list(object[["posteriors"]][["upsilon"]])) {
      if ("coeffs" %in% names(object[["posteriors"]][["upsilon"]])) {
        n_gamma_glo <- nrow(object[["posteriors"]][["upsilon"]][["coeffs"]])
      }
    } else {
      n_gamma_glo <- nrow(object[["posteriors"]][["gamma_glo"]])
    }
    if ((n_gamma_glo / tt) %% k^2 == 0 & n_gamma_glo / tt >= 1) {
      tvp_gamma_glo <- TRUE
      n_gamma_glo <- n_gamma_glo / tt
    }
    
    gamma_glo <- object[["posteriors"]][["upsilon"]]
    object[["posteriors"]][["upsilon"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(gamma_glo, tvp_gamma_glo, n_gamma_glo, tt, "gamma_global"))
    rm(gamma_glo)
  } else {
    object[["posteriors"]][["upsilon"]] <- NULL
  }
  
  # Gamma deterministic ----
  if(!is.null(object[["posteriors"]][["c"]])) {
    if (is.list(object[["posteriors"]][["c"]])) {
      if ("coeffs" %in% names(object[["posteriors"]][["c"]])) {
        n_gamma_det <- nrow(object[["posteriors"]][["c"]][["coeffs"]])
      }
    } else {
      n_gamma_det <- nrow(object[["posteriors"]][["c"]])
    }
    if ((n_gamma_det / tt) %% k^2 == 0 & n_gamma_det / tt >= 1) {
      tvp_gamma_det <- TRUE
      n_gamma_det <- n_gamma_det / tt
    }
    
    gamma_det <- object[["posteriors"]][["c"]]
    object[["posteriors"]][["c"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(gamma_det, tvp_gamma_det, n_gamma_det, tt, "gamma_deterministic"))
    rm(gamma_det)
  } else {
    object[["posteriors"]][["c"]] <- NULL
  }
  
  if(!is.null(object[["posteriors"]][["sigma"]])) {
    if (is.list(object[["posteriors"]][["sigma"]])) {
      if ("coeffs" %in% names(object[["posteriors"]][["sigma"]])) {
        n_sigma <- nrow(object[["posteriors"]][["sigma"]][["coeffs"]])
      }
    } else {
      n_sigma <- nrow(object[["posteriors"]][["sigma"]])
    }
    if ((n_sigma / tt) %% k == 0 & n_sigma / tt >= 1) {
      tvp_sigma <- TRUE
      n_sigma <- n_sigma / tt
    }
    if (n_sigma %% (k * k) != 0) {
      stop("Row number of coefficient draws of 'Sigma' is not k^2 or multiples thereof.")
    }
    
    sigma <- object[["posteriors"]][["sigma"]]
    object[["posteriors"]][["sigma"]] <- NULL
    object[["posteriors"]] <- c(object[["posteriors"]], .country_fill_helper(sigma, tvp_sigma, n_sigma, tt, "sigma"))
    rm(sigma)
  }
  
  vars <- c("a0",
            "alpha",
            "beta_domestic", "beta_foreign", "beta_global", "beta_determinisitc",
            "gamma_domestic", "gamma_foreign", "gamma_global", "gamma_determinisitc",
            "sigma")
  
  for (i in vars) {
    if (!is.null(object[["posteriors"]][[i]])) {
      attr(object[["posteriors"]][[i]], "mcpar") <- mc_spec
    }
  }
  
  class(object) <- append("ctryvecest", class(object))
  
  return(object)
}