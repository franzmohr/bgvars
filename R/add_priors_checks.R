
.add_priors_check_bvs <- function(object, bvs) {
  
  use_bvs_error <- FALSE
  
  if (object[["model"]][["varsel"]] == "bvs") {
    
    if (is.null(bvs)) {
      stop("BVS was chosen as variable selection algorithm, but no prior specification was provided.")
    }
    
    if (is.null(bvs[["inprior"]])) {
      stop("Argument 'bvs$inprior' must be specified for BVS")
    }
    
    if (!is.null(bvs[["covar"]])) {
      if (bvs[["covar"]]) {
        use_bvs_error <- TRUE 
      }
    }
  }
}

.add_priors_check_coef <- function(object, coef) {
  
  allowed_coef_arguments <- c("v_i", "v_i_det", "coint_var", "const", "minnesota",
                              "max_var", "shape", "rate", "rate_det", "omega_v")
  # Only an error correction sub-model has loadings to give a rate of their own.
  if (inherits(object, "vecxsubmodel")) {
    allowed_coef_arguments <- c(allowed_coef_arguments, "rate_alpha", "omega_v_alpha")
  }
  for (i in names(coef)) {
    if (!i %in% allowed_coef_arguments) {
      stop(paste0("Element '", i, "' in argument 'coef' is not recognised."))
    }
  }

  if (!is.null(coef[["v_i"]])) {
    if (coef[["v_i"]] < 0) {
      stop("Argument 'v_i' must be at least 0.")
    } 
    # Define "v_i_det" if not specified (needed for a check later)
    if (is.null(coef[["v_i_det"]])) {
      coef[["v_i_det"]] <- coef[["v_i"]]
    }
  } else {
    if (!any(c("minnesota", "ssvs") %in% names(coef))) {
      stop("If 'coef$v_i' is not specified, at least 'coef$minnesota' or 'coef$ssvs' must be specified.")
    }
  }
  
  # Tests for specifications used in TVP models
  if (!is.null(object[["model"]][["tvp"]])) {
    if (object[["model"]][["tvp"]]) {
      # The state variances take one of two parameterisations and never both.
      # "shape" and "rate" put a gamma prior on the precisions of the state
      # equation, the centred form; "omega_v" puts a normal prior on the signed
      # standard deviations instead, the non-centred form of
      # Fruehwirth-Schnatter and Wagner (2010). The samplers reject a block that
      # names both, so it is caught here where the argument is still visible.
      if ("omega_v" %in% names(coef)) {
        centred <- c("shape", "rate", "rate_det", "rate_alpha")
        given <- centred[centred %in% names(coef)]
        if (length(given) > 0) {
          stop("Argument 'coef' gives both a gamma prior on the precisions of the ",
               "state equation ('", paste(given, collapse = "', '"),
               "') and a normal prior on their standard deviations ('omega_v'). ",
               "Give one: 'omega_v' for the non-centred parameterisation, 'shape' ",
               "and 'rate' for the centred one.")
        }
        if (any(coef[["omega_v"]] <= 0)) {
          stop("Argument 'coef$omega_v' must be positive. It is the prior variance ",
               "of the standard deviation of the state innovations, and also the ",
               "prior mean of their variance.")
        }
        if (!is.null(coef[["omega_v_alpha"]]) && any(coef[["omega_v_alpha"]] <= 0)) {
          stop("Argument 'coef$omega_v_alpha' must be positive.")
        }
      } else {
        if ("omega_v_alpha" %in% names(coef)) {
          stop("Argument 'coef$omega_v_alpha' is the prior of the loadings under the ",
               "non-centred parameterisation of the state equation and needs ",
               "'coef$omega_v' beside it.")
        }
        if (!"shape" %in% names(coef)) {
          stop("Argument 'coef$shape' must be specified for TVP models, or ",
               "'coef$omega_v' for the non-centred parameterisation of the state ",
               "equation.")
        }
        if (!"rate" %in% names(coef)) {
          stop("Argument 'coef$rate' must be specified for TVP models, or ",
               "'coef$omega_v' for the non-centred parameterisation of the state ",
               "equation.")
        }
      }
      
      if (!is.null(coef[["const"]])) {
        if ("character" %in% class(coef[["const"]])) {
          if (!coef[["const"]] %in% c("first", "mean")) {
            stop("Invalid specificatin of coef$const.")
          }
        }
      }
    } 
  }
  
  if (!is.null(coef[["minnesota"]])) {
    if (!"list" %in% class(coef[["minnesota"]])) {
      stop("Argument coef$minnesota must be a named list.")
    }
    if (is.null(names(coef[["minnesota"]]))) {
      stop("Argument coef$minnesota must be a named list.")
    }
    if (!all(c("kappa1", "kappa2", "kappa4") %in% names(coef[["minnesota"]]))) {
      stop("Argument coeff$minnesota must contain at least the elements 'kappa1', 'kappa2' and 'kappa4'.")
    }
    if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & is.null(coef[["v_i"]])) {
      stop("If error covarances should be estimated, argument coef$v_i must be provided also when the Minnesota prior is used.")
    }
  }
}


.add_priors_check_sigma <- function(object, sigma) {
  
  if (is.null(sigma)) {
    stop("Argument 'sigma' may not be NULL.")
  }
  
  if (length(sigma) < 2) {
    stop("Argument 'sigma' must be at least of length 2.")
  } else {
    error_prior <- NULL
    
    if (object$model$error %in% c("gamma", "gamma+covar")) {
      if (all(c("shape", "rate") %in% names(sigma))) {
        error_prior <- "gamma"
      } else {
        stop("Gamma prior requires specification of elements 'shape' and 'rate' in 'sigma'.")
      }
      if (sigma$shape < 0) {
        stop("Argument 'sigma$shape' must be at least 0.")
      }
      if (sigma$rate <= 0) {
        stop("Argument 'sigma$rate' must be larger than 0.")
      } 
    }
    
    if (object$model$error %in% c("sv", "sv+covar")) {
      # The innovations of the log-volatilities take one of the two
      # parameterisations the coefficients take, and never both: "shape" and
      # "rate" put a gamma prior on their precision, "omega_v" a normal prior on
      # their signed standard deviation. The second is only read by a model whose
      # parameters vary over time, so it is refused elsewhere rather than written
      # into the priors and silently left unread.
      if (!is.null(sigma[["omega_v"]])) {
        if (!isTRUE(object[["model"]][["tvp"]])) {
          stop("Argument 'sigma$omega_v' is the non-centred parameterisation of the ",
               "state equation of the log-volatilities and is only available for ",
               "sub-models with time varying parameters. Use 'sigma$shape' and ",
               "'sigma$rate' for the centred one.")
        }
        given <- c("shape", "rate")[c("shape", "rate") %in% names(sigma)]
        if (length(given) > 0) {
          stop("Argument 'sigma' gives both a gamma prior on the precisions of the ",
               "log-volatility innovations ('", paste(given, collapse = "', '"),
               "') and a normal prior on their standard deviations ('omega_v'). ",
               "Give one: 'omega_v' for the non-centred parameterisation, 'shape' ",
               "and 'rate' for the centred one.")
        }
        if (!is.numeric(sigma[["omega_v"]]) || any(!is.finite(sigma[["omega_v"]])) ||
            any(sigma[["omega_v"]] <= 0)) {
          stop("Argument 'sigma$omega_v' must be a finite positive numeric. It is the ",
               "prior variance of the standard deviation of the log-volatility ",
               "innovations, and also the prior mean of their variance.")
        }
      }
      state_prior <- if (is.null(sigma[["omega_v"]])) c("shape", "rate") else "omega_v"
      if (any(!c("mu", "v_i", state_prior, "state_variance", "offset") %in% names(sigma))) {
        stop("Missing prior specifications for stochastic volatility prior.")
      }
      error_prior <- "sv"
    }
    
    if (object$model$error == "wishart") {
      if (all(c("df", "scale") %in% names(sigma))) {
        error_prior <- "wishart"
      } else {
        stop("Wishart prior requires specification of elements 'df' and 'scale' in 'sigma'.")
      }
      if (object[["model"]][["structural"]]) {
        stop("Structural models may not use a Wishart prior. Consider using a gamma prior instead.")
      }
      if (sigma$df < 0) {
        stop("Argument 'sigma$df' must be at least 0.")
      }
      if (sigma$scale <= 0) {
        stop("Argument 'sigma$scale' must be larger than 0.")
      } 
    }
    
    if (is.null(error_prior)) {
      stop("Invalid specification for argument 'sigma'.")
    }
  }
  
  return(error_prior)
}




.add_priors_check_ssvs <- function(object, ssvs) {
  
  use_ssvs_error <- FALSE
  
  if (is.null(ssvs[["inprior"]])) {
    stop("Argument 'varsel$inprior' must be specified for SSVS.")
  }
  if (is.null(ssvs[["tau"]]) & is.null(ssvs[["semiautomatic"]])) {
    stop("Either argument 'varsel$tau' or 'varsel$semiautomatic' must be specified for SSVS.")
  }
  if (!is.null(ssvs[["covar"]])) {
    if (ssvs[["covar"]]) {
      use_ssvs_error <- TRUE 
    }
  }
  if (object[["model"]][["error"]] == "gamma+covar" & use_ssvs_error & is.null(ssvs[["tau"]])) {
    stop("If SSVS should be applied to error covariances, argument 'varsel$tau' must be specified.")
  }
  if (object[["model"]][["structural"]] & is.null(ssvs[["tau"]])) {
    stop("If SSVS should be used with structural models, argument 'varsel$tau' must be specified.")
  }
  
  if (!is.null(ssvs[["semiautomatic"]])) {
    if (!"numeric" %in% class(ssvs[["semiautomatic"]])) {
      stop("Argument 'varsel$semiautomatic' must be a numeric vector.")
    }
    if (length(ssvs[["semiautomatic"]]) != 2) {
      stop("Argument 'varsel$semiautomatic' must be a numeric vector with two elements.")
    }
  }
}

# Warns when BVS is combined with a flat prior on the coefficients. A Minnesota
# prior is never flat, and it may come without v_i, which is why it is checked
# first. Without one, the precision of the deterministic terms is v_i_det, or
# v_i when that is not given -- the fallback the prior itself uses.
.add_priors_warn_bvs_flat <- function(coef, minnesota, exclude_det) {
  if (minnesota) {
    return(invisible(NULL))
  }
  v_i_det <- coef[["v_i_det"]]
  if (is.null(v_i_det)) {
    v_i_det <- coef[["v_i"]]
  }
  if (coef[["v_i"]] == 0 || (v_i_det == 0 && !exclude_det)) {
    warning("Using BVS with an uninformative prior is not recommended.")
  }
  invisible(NULL)
}
