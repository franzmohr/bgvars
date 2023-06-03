#' Add Priors to Country-Specific VARX Models of a GVAR Model
#'
#' Adds prior specifications to a list of country models, which was produced by
#' the function \code{\link{create_models}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_models}}.
#' @param coef a named list of prior specifications for the coefficients of the country-specific
#' models. For the default specification all prior means are set to zero and the diagonal elements of
#' the inverse prior variance-covariance matrix are set to 1 for coefficients corresponding to non-deterministic
#' terms. For deterministic coefficients the prior variances are set to 10 via \code{v_i_det = 0.1}.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' For further specifications see 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix
#' of the country models. For the default specification of an inverse Wishart distribution
#' the prior degrees of freedom are set to the number of endogenous variables and
#' the prior variances to 1. See 'Details'.
#' @param ssvs a named list of prior specifications for the SSVS algorithm. See 'Details'.
#' @param bvs a named list of prior specifications for the BVS algorithm. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The argument \code{coef} can contain the following elements
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients corresponding to deterministic terms. Default is 0.1.}
#'   \item{\code{coint_var}}{a logical specifying whether the prior mean of the first own lag of an
#'   endogenous variable in a VAR model should be set to 1. Default is \code{FALSE}.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{const = "mean"}, the means of the series of endogenous variables are used as prior means.
#'   If \code{const = "first"}, the first values of the series of endogenous variables are used as prior means.}
#'   \item{\code{minnesota}}{a list of length 4 containing parameters for the calculation of
#'   the Minnesota prior, where the element names are \code{kappa0}, \code{kappa1}, \code{kappa2} and \code{kappa3}.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th
#'   lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{0}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{0} \kappa_{1}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{0} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for foreign and global exogenous variables,}}
#'   \deqn{ \kappa_{0} \kappa_{3} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   If \code{kappa2 = NULL}, \eqn{\kappa_{0} \kappa_{3} \sigma_{i}^2} will be used for foreign
#'   and global exogenous variables instead.
#'   
#'   For VEC models the function only provides priors for the non-cointegration part of the model. The
#'   residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}}{a numeric specifying the maximum prior variance that is allowed for
#'   non-deterministic coefficients.}
#'   \item{\code{shape}}{an integer specifying the prior degrees of freedom of the error term of the state equation. Default is 3.}
#'   \item{\code{rate}}{a numeric specifying the prior error variance of the state equation. Default is 0.0001.}
#'   \item{\code{rate_det}}{a numeric specifying the prior error variance of the state equation corresponding to deterministic terms. Default is 0.0001.}
#' }
#' If \code{minnesota} is specified, \code{v_i} and \code{v_i_det} are ignored.
#' 
#' Argument \code{sigma} can contain the following elements:
#' \describe{
#'   \item{\code{df}}{an integer or character specifying the prior degrees of freedom of the error term. Only
#'   used, if the prior is inverse Wishart.
#'   Default is \code{"k"}, which indicates the amount of endogenous variables in the respective model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If an integer
#'   is provided, the degrees of freedom are set to this value in all models.}
#'   \item{\code{scale}}{a numeric specifying the prior error variance of endogenous variables.
#'   Default is 1.}
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error term. Only
#'   used, if the prior is inverse gamma or if time varying volatilities are estimated.
#'   For models with constant volatility the default is \code{"k"}, which indicates the amount of endogenous
#'   variables in the respective country model. \code{"k + 3"} can be used to set the prior to the amount of
#'   endogenous variables plus 3. If a numeric is provided, the shape parameters are set to this value in all
#'   models. For models with stochastic volatility this prior refers to the error variance of the state
#'   equation.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term. Only used, if the
#'   prior is inverse gamma or if time varying volatilities are estimated. For models with stochastic
#'   volatility this prior refers to the error variance of the state equation.}
#'   \item{\code{mu}}{numeric of the prior mean of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{v_i}}{numeric of the prior precision of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{sigma_h}}{numeric of the initial draw for the variance of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{covar}}{logical indicating whether error covariances should be estimated. Only used
#'   in combination with an inverse gamma prior or stochastic volatility, for which \code{shape} and
#'   \code{rate} must be specified.}
#' }
#' \code{df} and \code{scale} must be specified for an inverse Wishart prior. \code{shape} and \code{rate}
#' are required for an inverse gamma prior. For structural models or models with stochastic volatility
#' only a gamma prior specification is allowed.
#' 
#' The argument \code{ssvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model. Default is 0.5.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second. Default is \code{c(0.05, 10)}.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the country VARX model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables. This is the
#'   semiautomatic approach described in George et al. (2008).}
#'   \item{\code{exclude_det}}{logical indicating whether deterministic terms should be excepted from SSVS.}
#' }
#' Either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' The argument \code{bvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model. Default is 0.5.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters for the calculation of
#'   Minnesota-like inclusion priors, where
#'   \deqn{ \kappa_{1}{l} \textrm{ for own lags of domestic endogenous variables,}}
#'   \deqn{ \frac{\kappa_{2}}{l} \textrm{ for domestic endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{3}}{l + 1} \textrm{ for foreign and global exogenous variables,}}
#'   \deqn{ \kappa_{4} \textrm{ for deterministic terms.}}
#'   The indices of \eqn{\kappa} correspond to the positions of the elements in the argument \code{kappa}.}
#'   \item{\code{exclude_det}}{logical indicating whether deterministic terms should be excepted from BVS.}
#' }
#' 
#' @return A list of country models.
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
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230.
#' 
#' @export
add_priors.gvarsubmodels <- function(object,
                                     coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001, rate_det = 0.01),
                                     sigma = list(df = "k", scale = 1, mu = 0, v_i = 0.01, sigma_h = 0.05),
                                     ssvs = NULL,
                                     bvs = NULL,
                                     ...){
  
  # Checks - Coefficient priors ----
  if (!is.null(coef)) {
    if (!is.null(coef[["v_i"]])) {
      if (coef[["v_i"]] < 0) {
        stop("Argument 'v_i' must be at least 0.")
      } 
      if (is.null(coef[["v_i_det"]])) {
        coef[["v_i_det"]] <- coef[["v_i"]]
      }
    } else {
      if (!any(c("minnesota", "ssvs") %in% names(coef))) {
        stop("If 'coef$v_i' is not specified, at least 'coef$minnesota' or 'coef$ssvs' must be specified.")
      }
    }
  }
  
  if (!is.null(coef[["const"]])) {
    if (class(coef[["const"]]) == "character") {
      if (!coef[["const"]] %in% c("first", "mean")) {
        stop("Invalid specificatin of coef$const.")
      }
    }
  }
  
  # Checks - Error priors ----
  if (length(sigma) < 2) {
    stop("Argument 'sigma' must be at least of length 2.")
  } else {
    error_prior <- NULL
    if (any(unlist(lapply(object, function(x) {x[["model"]][["sv"]]})))) { # Check for SV
      if (any(!c("mu", "v_i", "shape", "rate") %in% names(sigma))) {
        stop("Missing prior specifications for stochastic volatility prior.")
      }
      error_prior <- "sv"
    } else {
      
      if (all(c("df", "scale", "shape", "rate") %in% names(sigma))) {
        error_prior <- "both"
      } else {
        if (all(c("shape", "rate") %in% names(sigma))) {
          error_prior <- "gamma"
        }
        if (all(c("df", "scale") %in% names(sigma))) {
          error_prior <- "wishart"
        } 
      }

      if (is.null(error_prior)) {
        stop("Invalid specification for argument 'sigma'.")
      }
      # if (error_prior == "wishart" & any(unlist(lapply(object, function(x) {x[["model"]][["structural"]]})))) {
      #   stop("Structural models may not use a Wishart prior. Consider specifying arguments 'sigma$shape' and 'sigma$rate' instead.")
      # }
      
      if (error_prior == "wishart" | error_prior == "both") {
        if (sigma[["df"]] < 0) {
          stop("Argument 'sigma$df' must be at least 0.")
        }
        if (sigma[["scale"]] <= 0) {
          stop("Argument 'sigma$scale' must be larger than 0.")
        } 
      }
      if (error_prior == "gamma" | error_prior == "both") {
        if (sigma[["shape"]] < 0) {
          stop("Argument 'sigma$shape' must be at least 0.")
        }
        if (sigma[["rate"]] <= 0) {
          stop("Argument 'sigma$rate' must be larger than 0.")
        } 
      }
    }
  }
  
  # Check Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  # Check coint VAR ----
  coint_var <- FALSE # Cointegrated VAR?
  if (!is.null(coef[["coint_var"]])) {
    if (coef[["coint_var"]]) {
      coint_var <- TRUE 
    }
  }
  
  # Check SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (!is.null(ssvs)) {
    
    if (error_prior == "sv") {
      stop("SSVS is not supported for stochastic volatility models. You could use BVS instead.")
    }
    
    if (is.null(ssvs[["inprior"]])) {
      stop("Argument 'ssvs$inprior' must be specified for SSVS.")
    }
    if (is.null(ssvs[["tau"]]) & is.null(ssvs[["semiautomatic"]])) {
      stop("Either argument 'ssvs$tau' or 'ssvs$semiautomatic' must be specified for SSVS.")
    }
    if (is.null(ssvs[["exclude_det"]])) {
      ssvs[["exclude_det"]] <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approch of 
    # George et al. (2008) should be used
    if (!is.null(ssvs[["semiautomatic"]])) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
    
    if (!is.null(ssvs[["covar"]])) {
      if (ssvs[["covar"]]) {
        if (error_prior == "wishart") {
          stop("If SSVS should be applied to error covariances, argument 'sigma$shape' must be specified.")
        }
        use_ssvs_error <- TRUE 
      }
      if (is.null(ssvs[["tau"]])) {
        stop("If SSVS should be applied to error covariances, argument 'ssvs$tau' must be specified.")
      }
    }
  }
  
  # Check BVS ---- 
  # Prior a la Korobilis 2013
  use_bvs <- FALSE
  use_bvs_error <- FALSE
  if (!is.null(bvs)) {
    use_bvs <- TRUE
    if (is.null(bvs[["inprior"]])) {
      stop("If BVS should be applied, argument 'bvs$inprior' must be specified.")
    }
    if (is.null(bvs[["exclude_det"]])) {
      bvs[["exclude_det"]] <- FALSE
    }
    if (!is.null(bvs[["covar"]])) {
      if (bvs[["covar"]]) {
        if (error_prior == "wishart") {
          stop("If BVS should be applied to error covariances, argument 'sigma$shape' must be specified.")
        }
        use_bvs_error <- TRUE 
      }
    }
    if (coef[["v_i"]] == 0 | (coef[["v_i_det"]] == 0 & !bvs[["exclude_det"]])) {
      warning("Using BVS with an uninformative prior is not recommended.")
    }
  }
  
  if (use_ssvs & use_bvs) {
    stop("SSVS and BVS cannot be applied at the same time.")
  }
  
  if (error_prior == "wishart" & (use_ssvs_error | use_bvs_error)) {
    stop("Wishart prior not allowed when BVS or SSVS are applied to covariance matrix.")
  }
  
  # Generate priors for each country ----
  for (i in 1:length(object)) {
    
    # Get model specs to obtain total number of coeffs
    k <- length(object[[i]][["model"]][["domestic"]][["variables"]])
    p_domestic <- object[[i]][["model"]][["domestic"]][["lags"]]
    
    if (k == 1 & (use_ssvs_error | use_bvs_error)) {
      stop("BVS or SSVS cannot be applied to covarianc matrix when there is only one endogenous variable.")
    }
    
    m <- length(object[[i]][["model"]][["foreign"]][["variables"]])
    p_foreign <- object[[i]][["model"]][["foreign"]][["lags"]]
    
    global <- !is.null(object[[i]][["model"]][["global"]])
    if (global) {
      n <- length(object[[i]][["model"]][["global"]][["variables"]])
      p_global <- object[[i]][["model"]][["global"]][["lags"]]
    } else {
      n <- 0
      p_global <- 0
    }
    
    # Add a lag to foreign and global model for VARs to include
    # contemporary variables
    p_foreign <- p_foreign + 1
    if (global) {
      p_global <- p_global + 1
    }
    
    # Total # of non-deterministic coefficients
    tot_par <- k * (k * p_domestic + m * p_foreign + n * p_global)
    
    covar <- FALSE
    if (!is.null(sigma[["covar"]])) {
      covar <- sigma[["covar"]]
    }
    structural <- object[[i]][["model"]][["structural"]]
    if (covar & structural) {
      stop("Error covariances and structural coefficients cannot be estimated at the same time.")
    }
    n_struct <- 0
    if (structural & k > 1) {
      n_struct <- (k - 1) * k / 2
      tot_par <- tot_par + n_struct
    }
    
    tvp <- object[[i]][["model"]][["tvp"]]
    sv <- object[[i]][["model"]][["sv"]]
    
    # Add number of non-cointegration deterministic terms
    n_det <- 0
    if (!is.null(object[[i]][["data"]][["deterministic"]])){
      n_det <- ncol(object[[i]][["data"]][["deterministic"]]) * k
      tot_par <- tot_par + n_det
    }
    
    # Priors ----
    ## Coefficients ----
    if (tot_par > 0) {
      
      #### Prior means ----
      
      mu <- matrix(rep(0, tot_par - n_struct), k)
      
      # Add 1 to first own lags for cointegrated VAR model
      if (coint_var & p_domestic > 0) {
        mu[1:k, 1:k] <- diag(1, k)
      }
      
      # Prior for intercept terms
      if (n_det > 0) {
        
        if (!is.null(coef[["const"]]))  {
          
          pos <- which(dimnames(object[[i]][["data"]][["Z"]])[[2]] == "const")
          
          if (length(pos) == 1) {
            if ("character" %in% class(coef[["const"]])) {
              if (coef[["const"]] == "first") {
                mu[, pos] <- object[[i]][["data"]][["Y"]][1, ]
              }
              if (coef[["const"]] == "mean") {
                mu[, pos] <- colMeans(object[[i]][["data"]][["Y"]])
              }
            }
            if ("numeric" %in% class(coef[["const"]])) {
              mu[, pos] <- coef[["const"]]
            } 
          }
        }
      }
      
      mu <- matrix(mu)
      
      if (structural) {
        mu <- rbind(mu, matrix(0, n_struct))
      }
      
      object[[i]][["priors"]][["coefficients"]] <- list(mu = mu)
      
      # Obtain OLS estimates for the calculation of the
      # Minnesota prior or the semi-automatic SSVS approach
      # and for use as initial values
      
      # Get data
      y <- t(object[[i]][["data"]][["Y"]])
      if (tot_par > 0 & tot_par < NCOL(y)) {
        z <- object[[i]][["data"]][["SUR"]]
        ols <- solve(crossprod(z)) %*% crossprod(z, matrix(y))
        u <- matrix(matrix(y) - z %*% ols, NROW(y))
        object[[i]][["initial"]][["coefficients"]][["draw"]] <- ols
      } else {
        if (tot_par > 0) {
          object[[i]][["initial"]][["coefficients"]][["draw"]] <- mu
        }
        u <- y - matrix(apply(y, 1, mean), nrow = NROW(y), ncol = NCOL(y))
      }
      x <- t(object[[i]][["data"]][["Z"]])
      tt <- ncol(y)
      ols_sigma <- y %*% (diag(1, tt) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / (tt - nrow(x))
      
      if (minnesota) {
        # Determine positions of deterministic terms for calculation of sigma
        pos_det <- NULL
        if (!is.null(object[[i]][["model"]][["deterministic"]])) {
          pos_det <- k * p_domestic + m * p_foreign + n * p_global + 1:length(object[[i]][["model"]][["deterministic"]])
        }
        
        # Obtain sigmas for V_i via estimation of AR model for each endogenous variable
        s_domestic <- diag(0, k)
        for (j in 1:k) {
          pos <- c(j + k * ((1:p_domestic) - 1), pos_det)
          y_temp <- matrix(y[j, ], 1)
          x_temp <- matrix(x[pos, ], length(pos))
          s_domestic[j, j] <- y_temp %*% (diag(1, tt) - t(x_temp) %*% solve(tcrossprod(x_temp)) %*% x_temp) %*% t(y_temp) / (tt - length(pos))
        }
        s_domestic <- sqrt(diag(s_domestic)) # Residual standard deviations (OLS)
        
        # Generate prior matrices ----
        
        # Minnesota prior ----
        V <- matrix(rep(NA, tot_par - n_struct), k) # Set up matrix for variances
        
        # Endogenous variables
        if (p_domestic > 0) {
          for (r in 1:p_domestic) {
            for (l in 1:k) {
              for (j in 1:k) {
                if (l == j) {
                  V[l, (r - 1) * k + j] <- coef[["minnesota"]][["kappa0"]] / r^2
                } else {
                  V[l, (r - 1) * k + j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa1"]] / r^2 * s_domestic[l]^2 / s_domestic[j]^2
                }
              }
            }
          }
        }
        
        # Weakly exogenous variables
        s_foreign <- sqrt(apply(matrix(x[k * p_domestic + 1:m, ], m), 1, stats::var))
        for (r in 1:p_foreign) {
          for (l in 1:k) {
            for (j in 1:m) {
              # Note that r starts at 1, so that this is equivalent to l + 1
              if (is.null(coef[["minnesota"]][["kappa2"]])) {
                V[l, p_domestic * k + (r - 1) * m + j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa3"]] * s_domestic[l]^2
              } else {
                V[l, p_domestic * k + (r - 1) * m + j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa2"]] / r^2 * s_domestic[l]^2 / s_foreign[j]^2
              }
            }
          }
        }
        
        # Glogal variables
        if (global) {
          s_global <- sqrt(apply(matrix(x[k * p_domestic + m * p_foreign + 1:n, ], n), 1, stats::var))
          for (r in 1:p_global) {
            for (l in 1:k) {
              for (j in 1:n) {
                # Note that r starts at 1, so that this is equivalent to l + 1
                if (is.null(coef[["minnesota"]][["kappa2"]])) {
                  V[l, p_domestic * k + m * p_foreign + (r - 1) * n + j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa3"]] * s_domestic[l]^2
                } else {
                  V[l, p_domestic * k + m * p_foreign + (r - 1) * n + j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa2"]] / r^2 * s_domestic[l]^2 / s_global[j]^2
                }
              }
            }
          }
        }
        
        # Restrict prior variances
        if (!is.null(coef[["max_var"]])) {
          if (any(stats::na.omit(V) > coef[["max_var"]])) {
            V[which(V > coef[["max_var"]])] <- coef[["max_var"]]
          }
        }
        
        # Deterministic variables
        if (!is.null(object[[i]][["data"]][["deterministic"]])){
          V[, -(1:(k * p_domestic + m * p_foreign + n * p_global))] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa3"]] * s_domestic^2
        }
        
        V <- matrix(V)
        
        # Structural parameters
        if (structural) {
          V_struct <- matrix(NA, k, k)
          for (j in 1:(k - 1)) {
            V_struct[(j + 1):k, j] <- coef[["minnesota"]][["kappa0"]] * coef[["minnesota"]][["kappa1"]] * s_domestic[(j + 1):k]^2 / s_domestic[j]^2  
          }
          V_struct <- matrix(V_struct[lower.tri(V_struct)])
          V <- rbind(V, V_struct)
        }
        
        v_i <- diag(c(1 / V))
        
        object[[i]][["priors"]][["coefficients"]][["v_i"]] <- v_i
      } # End of minnesota condition
      
      
      # Inclusion priors
      if (use_ssvs | use_bvs) {
        inprior <- matrix(NA, k, (tot_par - n_struct) / k)
        include <- matrix(1:tot_par)
        if (use_ssvs) {
          prob <- ssvs[["inprior"]]
          kappa <- ssvs[["minnesota"]]
          exclude_det <- ssvs[["exclude_det"]]
        }
        if (use_bvs) {
          prob <- bvs[["inprior"]]
          kappa <- bvs[["minnesota"]]
          exclude_det <- bvs[["exclude_det"]]
        }
        # For Minnesota-like inclusion parameters
        if (!is.null(kappa)) {
          # Domestic
          if (p_domestic > 0) {
            for (r in 1:p_domestic) {
              inprior[, (r - 1) * k + 1:k] <- kappa[2] / r
              if (k > 1) {
                diag(inprior[, (r - 1) * k + 1:k]) <- kappa[1] / r
              } else {
                inprior[, (r - 1) * k + 1] <- kappa[1] / r
              }
            }
          }
          
          # Foreign
          if (m > 0) {
            inprior[, p_domestic * k + 1:m] <- kappa[3]
            if (p_foreign > 0) {
              for (r in 1:p_foreign) {
                inprior[, p_domestic * k + (r - 1) * m + 1:m] <- kappa[3] / (1 + r)
              }
            }
          }
          
          # Global
          if (global) {
            inprior[, p_domestic * k + p_foreign * m + 1:n] <- kappa[3]
            if (p_global > 0) {
              for (r in 1:p_global) {
                inprior[, p_domestic * k + p_foreign * m + (r - 1) * n + 1:n] <- kappa[3] / (1 + r)
              }
            }
          }
          
          if (n_det > 0) {
            inprior[, p_domestic * k + p_foreign * m + p_global * n + 1:(n_det / k)] <- kappa[4]
          }
        } else {
          inprior[,] <- prob
        }
        
        inprior <- matrix(inprior)
        if (structural & k > 1) {
          inprior <- rbind(inprior, matrix(prob, n_struct))
        }
        
        if (n_det > 0 & exclude_det) {
          pos_det <- tot_par - n_det - n_struct + 1:n_det
          include <- matrix(include[-pos_det])
        }
      }
      
      # SSVS
      if (use_ssvs) {
        if (use_ssvs_semi) {
          # Semiautomatic approach
          cov_ols <- kronecker(solve(tcrossprod(x)), ols_sigma)
          se_ols <- sqrt(diag(cov_ols)) # OLS standard errors
          se_ols <- matrix(se_ols)
          
          tau0 <- se_ols * ssvs[["semiautomatic"]][1] # Prior if excluded
          tau1 <- se_ols * ssvs[["semiautomatic"]][2] # Prior if included
          
          if (structural) {
            warning("Semiautomatic approach for SSVS not available for structural variables. Using values of argument 'ssvs$tau' instead.")
            tau0 <- rbind(tau0, matrix(ssvs$tau[1], n_struct))
            tau1 <- rbind(tau1, matrix(ssvs$tau[2], n_struct))
          }
        } else {
          tau0 <- matrix(rep(ssvs[["tau"]][1], tot_par))
          tau1 <- matrix(rep(ssvs[["tau"]][2], tot_par))
        }
        
        object[[i]][["model"]][["varselect"]] <- "SSVS"
        
        object[[i]][["priors"]][["coefficients"]][["v_i"]] <- diag(1 / tau1[, 1]^2, tot_par)
        object[[i]][["priors"]][["coefficients"]][["ssvs"]][["inprior"]] <- inprior
        object[[i]][["priors"]][["coefficients"]][["ssvs"]][["include"]] <- include
        object[[i]][["priors"]][["coefficients"]][["ssvs"]][["tau0"]] <- tau0
        object[[i]][["priors"]][["coefficients"]][["ssvs"]][["tau1"]] <- tau1
      }
    }
    
    # Regular prior ----
    if (!minnesota & !use_ssvs) {
      v_i <- diag(coef[["v_i"]], tot_par)
      # Add priors for deterministic terms if they were specified
      if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
        diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
      }
      object[[i]][["priors"]][["coefficients"]][["v_i"]] <- v_i
    }
    
    # Prior variances of the state equations
    if (tvp) {
      object[[i]][["priors"]][["coefficients"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
      object[[i]][["priors"]][["coefficients"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
      if (n_det > 0 & !is.null(coef[["rate_det"]])) {
        object[[i]][["priors"]][["coefficients"]][["rate"]][tot_par - n_struct - n_det + 1:n_det] <- coef[["rate_det"]]
      }
    }
    
    # BVS
    if (use_bvs) {
      object[[i]][["model"]][["varselect"]] <- "BVS"
      object[[i]][["priors"]][["coefficients"]][["bvs"]][["inprior"]] <- inprior
      object[[i]][["priors"]][["coefficients"]][["bvs"]][["include"]] <- include
    }
    
    ## Covar priors ----
    
    if (!structural & covar & k > 1) {
      
      if (is.null(coef[["v_i"]])) {
        stop("If error covariances should be estimated, argument 'coef$v_i' must be specified.")
      }
      
      n_covar <- k * (k - 1) / 2
      object[[i]][["priors"]][["psi"]][["mu"]] <- matrix(0, n_covar)
      object[[i]][["priors"]][["psi"]][["v_i"]] <- diag(coef[["v_i"]], n_covar)
      if (object[[i]][["model"]][["tvp"]]) {
        object[[i]][["priors"]][["psi"]][["shape"]] <- matrix(coef[["shape"]], n_covar)
        object[[i]][["priors"]][["psi"]][["rate"]] <- matrix(coef[["rate"]], n_covar) 
      }
      
      # SSVS priors
      if (use_ssvs_error) {
        object[[i]][["priors"]][["psi"]][["ssvs"]][["inprior"]] <- matrix(ssvs[["inprior"]], n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["include"]] <- matrix(1:n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["tau0"]] <- matrix(ssvs[["tau"]][1], n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["tau1"]] <- matrix(ssvs[["tau"]][2], n_covar)
      }
      
      # BVS priors
      if (use_bvs_error) {
        object[[i]][["priors"]][["psi"]][["bvs"]][["inprior"]] <- matrix(bvs[["inprior"]], n_covar)
        object[[i]][["priors"]][["psi"]][["bvs"]][["include"]] <- matrix(1:n_covar)
      }
    }
    
    ## Error term ----
    if (object[[i]][["model"]][["sv"]]) {
      
      object[[i]][["priors"]][["sigma"]][["mu"]] <- matrix(sigma[["mu"]], k)
      object[[i]][["priors"]][["sigma"]][["v_i"]] <- diag(sigma[["v_i"]], k)
      object[[i]][["priors"]][["sigma"]][["shape"]] <- matrix(sigma[["shape"]], k)
      object[[i]][["priors"]][["sigma"]][["rate"]] <- matrix(sigma[["rate"]], k)
      
    } else {
      
      if (error_prior == "wishart" | (error_prior == "both" & !structural)) {
        object[[i]][["priors"]][["sigma"]][["type"]] <- "wishart"
        help_df <- sigma[["df"]]
        object[[i]][["priors"]][["sigma"]][["df"]] <- NA_real_
        object[[i]][["priors"]][["sigma"]][["scale"]] = diag(sigma[["scale"]], k)
      }
      
      if (error_prior == "gamma" | (error_prior == "both" & structural)) {
        object[[i]][["priors"]][["sigma"]][["type"]] <- "gamma"
        help_df <- sigma[["shape"]]
        object[[i]][["priors"]][["sigma"]][["shape"]] <- NA_real_
        object[[i]][["priors"]][["sigma"]][["rate"]] = matrix(sigma[["rate"]], k)
      }
      
      if (minnesota) {
        # Store LS estimate of variance coviariance matrix for analytical solution
        object[[i]][["priors"]][["sigma"]][["sigma_i"]] = solve(ols_sigma)
      }
      
      if (class(help_df) == "character") {
        if (grepl("k", help_df)) {
          # Transform character specification to expression and evaluate
          help_df <- eval(parse(text = help_df))
        } else {
          stop("Use no other letter than 'k' in 'sigma$df' to indicate the number of endogenous variables.")
        }
      }
      
      if (help_df < 0) {
        stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
      }
      
      if (error_prior == "wishart" | (error_prior == "both" & !structural)) {
        object[[i]][["priors"]][["sigma"]][["df"]] <- help_df
      }
      if (error_prior == "gamma" | (error_prior == "both" & structural)) {
        object[[i]][["priors"]][["sigma"]][["shape"]] <- matrix(help_df, k)
      }
    }
    
    # Initial values ----
    if (object[[i]][["model"]][["tvp"]]) {
      object[[i]][["initial"]][["coefficients"]][["sigma_i"]] <- diag(c(1 / object[[i]][["priors"]][["coefficients"]][["rate"]]), tot_par)
    }
    if (covar) {
      y_covar <- kronecker(-t(u), diag(1, k))
      pos <- NULL
      for (j in 1:k) {pos <- c(pos, (j - 1) * k + 1:j)}
      y_covar <- y_covar[, -pos]
      psi <- solve(crossprod(y_covar)) %*% crossprod(y_covar, matrix(u))
      object[[i]][["initial"]][["psi"]][["draw"]] <- psi
      Psi <- diag(1, k)
      for (j in 2:k) {
        Psi[j, 1:(j - 1)] <- t(psi[((j - 2) * (j - 1) / 2) + 1:(j - 1), 1])
      }
      u <- Psi %*% u
    }
    u <- apply(u, 1, stats::var)
    if (object[[i]][["model"]][["sv"]]) {
      object[[i]][["initial"]][["sigma"]][["h"]] <- log(matrix(u, nrow = NCOL(y), ncol = NROW(y), byrow = TRUE))
      object[[i]][["initial"]][["sigma"]][["sigma_h"]] <- matrix(sigma[["sigma_h"]], NROW(y))
    } else {
      object[[i]][["initial"]][["sigma"]][["sigma_i"]] <- diag(1 / u, NROW(y))
      dimnames(object[[i]][["initial"]][["sigma"]][["sigma_i"]]) <- list(dimnames(y)[[1]], dimnames(y)[[1]])
    }
  }
  return(object)
}
