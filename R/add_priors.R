#' Add Priors to Country Models
#'
#' Adds prior specifications to a list of country models, which was produced by
#' the function \code{\link{create_models}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_models}}.
#' @param coef a named list of prior specifications for the coefficients of the country-specific
#' models. For the default specification all prior means are set to zero and the diagonal elements of
#' the prior variance-covariance matrix are set to 1 for coefficients corresponding to non-deterministic
#' terms. For deterministic coefficients the prior variances are set to 10. The variances need to be
#' specified as precisions, i.e. as inverses of the variances. For further specifications see 'Details'.
#' @param coint a named list of prior specifications for cointegration coefficients of the
#' country-specific models. See 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix
#' of the country-specific models. For the default specification of an inverted Wishart distribution
#' the prior degrees of freedom are set to 3 and the prior variances to 1. See 'Details'.
#' @param ssvs a named list of prior specifications for the SSVS algorithm. See 'Details'.
#' @param bvs a named list of prior specifications for the BVS algorithm. See 'Details'.
#' 
#' @details The argument \code{coef} can contain the following elements
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients corresponding to deterministic terms. Default is 0.1.}
#'   \item{\code{coint_var}}{a logical specifying whether the prior mean of the first own lag of an
#'   endogenous variable should be set to 1. Default is \code{FALSE}.}
#'   \item{\code{kappa}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota prior. For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th
#'   lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of domestic endogenous variables,}}
#'   \deqn{ \frac{\kappa_{2}}{l^2} \frac{\sigma^2_{i}}{\sigma^2_{j}} \textrm{ for domestic endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{3}}{(l + 1)^2} \frac{\sigma^2_{i}}{\sigma^2_{j}} \textrm{ for foreign and global exogenous variables,}}
#'   \deqn{ \kappa_{4} \sigma^2_{i}  \textrm{ for deterministic terms.}}
#'   The indices of \eqn{\kappa} correspond to the positions of the elements in the argument \code{kappa}.
#'   \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   OLS estimate of the model. For exogenous variables \eqn{\sigma_{i}} corresponds to the standard
#'   deviation of the original series.
#'   
#'   For VEC models the function only provides priors for the non-cointegration part of the model. The
#'   residual standard errors \eqn{\sigma_i} are based on an unrestricted OLS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}}{numeric specifying the maximum prior variance that is allowed for non-deterministic coefficients.}
#' }
#' If \code{kappa} is specified, \code{v_i} and \code{v_i_det} will be ignored.
#' 
#' The argument \code{coint} can contain the following elements
#' \describe{
#'   \item{\code{coint_v_i}}{numeric between 0 and 1 specifying the shrinkage of the cointegration space prior. Default is 0.}
#'   \item{\code{coint_p_tau_i}}{numeric of the diagonal elements of the inverted prior matrix of
#' the central location of the cointegration space \eqn{sp(\beta)}. Default is 1.}
#' }
#' 
#' The argument \code{sigma} can contain the following elements
#' \describe{
#'   \item{\code{df}}{an integer specifying the prior degrees of freedom of the error term.
#'   Default is 3. If a VEC model is estimated, the rannk \eqn{r} of the cointegration matrix
#'   is automatically added to \code{df}.}
#'   \item{\code{scale}}{a numeric specifying the prior error variance of the endogenous variables.
#'   Default is 1.}
#' }
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
#' }
#' Either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' The argument \code{bvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model. Default is 0.5.}
#'   \item{\code{kappa}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota-like inclusion priors, where
#'   \deqn{ \kappa_{1} \textrm{ for own lags of domestic endogenous variables,}}
#'   \deqn{ \frac{\kappa_{2}}{l} \textrm{ for domestic endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{3}}{l + 1} \textrm{ for foreign and global exogenous variables,}}
#'   \deqn{ \kappa_{4} \textrm{ for deterministic terms.}}
#'   The indices of \eqn{\kappa} correspond to the positions of the elements in the argument \code{kappa}.}
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
#' # Generate EA area region with 2 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year rolling window averages
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
#' @export
add_priors <- function(object,
                       coef = list(v_i = 1, v_i_det = 0.1),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 1),
                       ssvs = list(inprior = 0.5, tau = c(0.05, 10)),
                       bvs = list(inprior = 0.5)){
  
  if (coef$v_i < 0) {
    stop("Argument 'v_i' must be at least 0.")
  }
  
  if (sigma$df < 0) {
    stop("Argument 'df' must be at least 0.")
  }
  
  if (sigma$scale < 0) {
    stop("Argument 'scale' must be at least 0.")
  }
  
  if (!is.null(bvs$kappa)) {
    if (length(bvs$kappa) != 4) {
      stop("Argument 'bvs$kappa' must have length 4.")
    }
    
    if (any(bvs$kappa > 1)) {
      stop("No element of argument 'bvs$kappa' may be larger than 1.")
    }
  }
  
  minnesota <- FALSE
  if (!is.null(coef$kappa)) {
    minnesota <- TRUE
  }
  
  coint_var <- FALSE
  if (!is.null(coef$coint_var)) {
    coint_var <- TRUE
  }
  
  for (i in 1:length(object)) {
    k <- length(object[[i]]$model$domestic$variables)
    p_domestic <- object[[i]]$model$domestic$lags
    m <- length(object[[i]]$model$foreign$variables)
    p_foreign <- object[[i]]$model$foreign$lags
    
    global <- !is.null(object[[i]]$model$global)
    if (global) {
      n <- length(object[[i]]$model$global$variables)
      p_global <- object[[i]]$model$global$lags
    } else {
      n <- 0
      p_global <- 0
    }
    
    if (object[[i]]$model$type == "VEC") {
      p_domestic <- p_domestic - 1
    } else {
      p_foreign <- p_foreign + 1
      if (global) {
        p_global <- p_global + 1
      }
    }
    
    # Total coefficients
    tot_par <- k * (k * p_domestic + m * p_foreign + n * p_global)
    
    # Add number of non-cointegration deterministic terms
    n_det <- 0
    if (!is.null(object[[i]]$data$deterministics)){
      if (object[[i]]$model$type == "VAR") {
        n_det <- ncol(object[[i]]$data$deterministics) * k
      }
      if (object[[i]]$model$type == "VEC") {
        if (!is.null(object[[i]]$data$deterministics$unrestricted)) {
          n_det <- ncol(object[[i]]$data$deterministics$unrestricted) * k
        }
      }
      tot_par <- tot_par + n_det
    }
    
    #### Cointegration (constant) ####
    if (object[[i]]$model$type == "VEC") {
      n_ect <- k * (k + m)
      if (global) {
        n_ect <- n_ect + k * n
      }
      if (!is.null(object[[i]]$data$deterministics$restricted)) {
        n_ect <- n_ect + k * ncol(object[[i]]$data$deterministics$restricted)
      }

      r_temp <- object[[i]]$model$cointegration$rank
      if (r_temp > 0) {
        object[[i]]$priors$cointegration <- list(v_i = coint$v_i,
                                                 p_tau_i = diag(coint$p_tau_i, n_ect / k))
      }
    }
    
    #### Non-cointegration (constant) ####
    ssvs_semi <- FALSE
    if (object[[i]]$model$variable_selection$type == "ssvs") {
      if (!is.null(ssvs$semiautomatic)) {
        ssvs_semi <- TRUE
      }
    }

    if (minnesota | ssvs_semi) {
      # Get data
      if (object[[i]]$model$type == "VAR") {
        temp <- .gen_varx(object[[i]])
        y <- temp$y
        x <- temp$x
      }
      if (object[[i]]$model$type == "VEC") {
        temp <- .gen_vecx(object[[i]])
        y <- temp$y
        x <- rbind(temp$ect, temp$x)
      }
      
      var_ols <- y %*% (diag(1, NCOL(y)) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / (NCOL(y) - nrow(x))
      s_domestic <- sqrt(diag(var_ols)) # Residual standard deviations (OLS)
    }
    
    if (tot_par > 0) {
      # Prior means
      mu <- matrix(rep(0, tot_par), k)
      if (coint_var & object[[i]]$model$type == "VAR") {
        mu[1:k, 1:k] <- diag(1, k)
      }
      mu <- matrix(mu)
      object[[i]]$priors$coefficients <- list(mu = mu)
      
      # Covariances
      if (minnesota) {

        V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
        
        # Endogenous variables
        if (p_domestic > 0) {
          for (r in 1:p_domestic) {
            for (l in 1:k) {
              for (j in 1:k) {
                if (l == j) {
                  V[l, (r - 1) * k + j] <- kappa[1] / r^2
                } else {
                  V[l, (r - 1) * k + j] <- kappa[2] / r^2 * s_domestic[l]^2 / s_domestic[j]^2
                }
              } 
            }
          } 
        }
        
        # Weakly exogenous variables
        if (object[[i]]$model$type == "VAR") {
          x <- t(x[k * p_domestic + 1:m, ])
        } else {
          x <- t(x[n_ect / k + k * p_domestic + 1:m, ])
        }
        s_foreign <- sqrt(apply(x, 2, stats::var))
        for (r in 1:p_foreign) {
          for (l in 1:k) {
            for (j in 1:m) {
              # Note that r starts at 1, so that this is equivalent to l + 1
              V[l, p_domestic * k + (r - 1) * m + j] <- kappa[3] / r^2 * s_domestic[l]^2 / s_foreign[j]^2
            }
          }
        }
        
        # Glogal variables
        if (global) {
          if (object[[i]]$model$type == "VAR") {
            x <- t(x[k * p_domestic + m * p_foreign + 1:n, ])
          } else {
            x <- t(x[n_ect / k + k * p_domestic + m * p_foreign + 1:n, ])
          }
          s_global <- sqrt(apply(x, 2, stats::var))
          rm(x)
          for (r in 1:p_global) {
            for (l in 1:k) {
              for (j in 1:n) {
                # Note that r starts at 1, so that this is equivalent to l + 1
                V[l, p_domestic * k + m * p_foreign + (r - 1) * n + j] <- kappa[3] / r^2 * s_domestic[l]^2 / s_foreign[j]^2
              }
            }
          } 
        }
        
        # Restrict prior variances
        if (!is.null(coef$max_var)) {
          if (any(stats::na.omit(V) > coef$max_var)) {
            V[which(V > coef$max_var)] <- coef$max_var
          } 
        }
        
        # Deterministic variables
        if (!is.null(object[[i]]$data$deterministics)){
          V[, -(1:(k * p_domestic + m * p_foreign + n * p_global))] <- kappa[4] * s_domestic^2
        }
        
        # Prior precision
        V_i_x <- diag(c(1 / V))
        object[[i]]$priors$coefficients$v_i <- V_i_x
      } else {
        v_i <- diag(coef$v_i, tot_par)
        if (n_det > 0) {
          diag(v_i)[tot_par - n_det + 1:n_det] <- coef$v_i_det
        }
        object[[i]]$priors$coefficients$v_i <- v_i
      }
    }
    
    #### Error term ####
    if (object[[i]]$model$type == "VEC") {
      if (!is.na(object[[i]]$model$cointegration$rank)) {
        df_temp <- sigma$df + r_temp 
      }
    } else {
      df_temp <- sigma$df
    }
    object[[i]]$priors$sigma <- list(df = df_temp,
                                     scale = diag(sigma$scale, k))
    
    #### Variable selection ####
    if (object[[i]]$model$variable_selection$type != "none") {
      # SSVS
      if (object[[i]]$model$variable_selection$type == "ssvs") {
        if (ssvs_semi) {
          # Semiautomatic approach
          cov_ols <- kronecker(solve(tcrossprod(x)), var_ols)
          se_ols <- matrix(sqrt(diag(cov_ols))) # OLS standard errors
          
          tau0 <- se_ols * ssvs$semiautomatic[1] # Prior if excluded
          tau1 <- se_ols * ssvs$semiautomatic[2] # Prior if included
        } else {
          tau0 <- matrix(rep(ssvs$tau[1], tot_par))
          tau1 <- matrix(rep(ssvs$tau[2], tot_par))
        }
        
        object[[i]]$priors$variable_selection$inprior <- matrix(ssvs$inprior, tot_par)
        object[[i]]$priors$variable_selection$tau0 <- tau0
        object[[i]]$priors$variable_selection$tau1 <- tau1
      }
      
      if (object[[i]]$model$variable_selection$type == "bvs") {
        if (is.null(bvs$kappa)) {
          object[[i]]$priors$variable_selection$inprior <- matrix(bvs$inprior, tot_par)
        } else {
          V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
          
          # Endogenous variables
          if (p_domestic > 0) {
            for (r in 1:p_domestic) {
              for (l in 1:k) {
                for (j in 1:k) {
                  if (l == j) {
                    V[l, (r - 1) * k + j] <- bvs$kappa[1]
                  } else {
                    V[l, (r - 1) * k + j] <- bvs$kappa[2] / r
                  }
                } 
              }
            } 
          }
          
          # Weakly exogenous variables
          for (r in 1:p_foreign) {
            for (l in 1:k) {
              for (j in 1:m) {
                # Note that r starts at 1, so that this is equivalent to l + 1
                V[l, p_domestic * k + (r - 1) * m + j] <- bvs$kappa[3] / r
              }
            }
          }
          
          # Glogal variables
          if (global) {
            for (r in 1:p_global) {
              for (l in 1:k) {
                for (j in 1:n) {
                  # Note that r starts at 1, so that this is equivalent to l + 1
                  V[l, p_domestic * k + m * p_foreign + (r - 1) * n + j] <- bvs$kappa[3] / r
                }
              }
            } 
          }
          
          # Deterministic variables
          if (n_det > 0){
            V[, -(1:(k * p_domestic + m * p_foreign + n * p_global))] <- bvs$kappa[4]
          }
          
          object[[i]]$priors$variable_selection$inprior <- matrix(V, tot_par) 
        }
      }
      
      if (object[[i]]$model$variable_selection$exclude_deterministics) {
        tot_par <- tot_par - n_det
      }
      object[[i]]$model$variable_selection$include <- matrix(1:tot_par, tot_par) 
    }
  }
  
  return(object)
}