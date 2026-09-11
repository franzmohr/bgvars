#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a VECX sub-model.
#' 
#' @param object an object of class 'vecxsubmodel', usually, a result of a call
#' to \code{\link{create_vecxsubmodel}}.
#' @param kappa1 a numeric specifying the prior variance of coefficients that correspond to
#' own lags of endogenous variables. See 'Details'.
#' @param kappa2 a numeric specifying the size of the prior variance of endogenous
#' variables, which do not correspond to own lags. See 'Details'.
#' @param kappa3 a numeric specifying the size of the prior variance of non-deterministic exogenous
#' variables. Default is \code{NULL}, which indicates that the formula
#' for the calculation of the prior variance of deterministic terms is used
#' for all exogenous variables. See 'Details'.
#' @param kappa4 a numeric specifying the size of the prior variance of deterministic
#' terms. See 'Details'.
#' @param max_var a positive numeric specifying the maximum prior variance that is allowed for
#' coefficients of non-deterministic variables. If \code{NULL} (default), the prior variances are not limited.
#' @param sigma either \code{"AR"} (default) or \code{"VAR"} indicating that the variances of the endogenous
#' variables \eqn{\sigma^2} are calculated based on a univariate AR regression or a least squares estimate of
#' the VAR form, respectively. In both cases all deterministic variables are used in the regressions,
#' if they appear in the model.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The function calculates the Minnesota prior in the same manner as for
#' a traditional VAR model. For the endogenous variable
#' \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#' \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#' \deqn{ \frac{\kappa_{1} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#' \deqn{ \frac{\kappa_{1} \kappa_{3}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for weakly exogenous and global variables,}}
#' \deqn{ \kappa_{1} \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#' where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#' LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#' In case structural parameters are estimated, the formula
#' \eqn{\kappa_{1} \kappa_{2} \frac{\sigma_{i}^2}{\sigma_{j}^2}} is used.
#' If \eqn{kappa_{3}} is not provided, prior variances are calculated in the same way as for
#' deterministic terms.
#' If the model does not contain exogenous variables, argument \code{kappa3} will be ignored.
#' 
#' The function only provides priors for the non-cointegration part of the model. The
#' residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#' endogenous variables on the error correction term and the non-cointegration regressors.
#' 
#' @return A list containing a matrix of prior means and the precision matrix of the cofficients and the
#' inverse variance-covariance matrix of the error term, which was obtained by an LS estimation.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2020). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
minnesota_prior.vecxsubmodel <- function(object, kappa1 = 2, kappa2 = .5, kappa3 = NULL, kappa4 = 5,
                                         max_var = NULL, sigma = "AR", ...) {
  
  if (kappa1 <= 0) {
    stop("Argument 'kappa1' must be positive.")
  }
  if (kappa2 <= 0) {
    stop("Argument 'kappa2' must be positive.")
  }
  if (!is.null(kappa3)) {
    if (kappa3 <= 0) {
      stop("Argument 'kappa3' must be positive.")
    } 
  }
  if (kappa4 <= 0) {
    stop("Argument 'kappa4' must be positive.")
  }
  if (!is.null(max_var)) {
    if (max_var <= 0) {
      stop("Argument 'max_var' must be positive.")
    } 
  }
  if (!sigma %in% c("AR", "VAR")) {
    stop("Argument 'sigma' must be either 'AR' or 'VAR'.")
  }
  
  y <- t(object[["data"]][["train"]][["y"]])
  k <- object[["model"]][["k"]]

  mu <- NULL  
  V <- NULL
  result <- NULL
  
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    
    if (!is.null(object[["data"]][["train"]][["x"]])) {
      
      
      x <- t(cbind(object[["data"]][["train"]][["w"]], object[["data"]][["train"]][["x"]]))
      n_ect <- NCOL(object[["data"]][["train"]][["w"]])
      tt <- NCOL(y)
      tot_par <- k * NROW(x)
      p_endogen <- object[["model"]][["p_endogen"]]
      p_endogen <- p_endogen - 1
      k_exogen <- object[["model"]][["k_exogen"]]
      p_exogen <- object[["model"]][["p_exogen"]]
      m <- object[["model"]][["m_global"]]
      s <- object[["model"]][["s_global"]]
      
      V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
      
      # Obtain OLS sigma
      ols_sigma <- y %*% (diag(1, tt) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / (tt - nrow(x))
      
      # Determine positions of deterministic terms for calculation of sigma
      pos_det <- NULL
      if (object[["model"]][["n_restricted"]] > 0 | object[["model"]][["n"]] > 0) {
        if (object[["model"]][["n_restricted"]] > 0) {
          pos_det <- c(pos_det, k + k_exogen + m + 1:length(object[["model"]][["n_restricted"]]))
        }
        if (object[["model"]][["n"]] > 0) {
          pos_det <- c(pos_det, n_ect + k * p_endogen + k_exogen * p_exogen + m * s + 1:length(object[["model"]][["n"]]))
        }
      }
      
      # Obtain sigmas for V_i
      if (sigma == "AR") { # Univariate AR
        s_endo <- diag(0, k)
        if (p_endogen > 0 | !is.null(pos_det)) {
          for (i in 1:k) {
            
            if (p_endogen > 0) {
              pos <- c(i, n_ect + i + k * ((1:p_endogen) - 1), pos_det)
            } else {
              pos <- c(i, pos_det)
            }
            
            y_temp <- matrix(y[i, ], 1)
            x_temp <- matrix(x[pos,], length(pos))
            s_endo[i, i] <- y_temp %*% (diag(1, tt) - t(x_temp) %*% solve(tcrossprod(x_temp)) %*% x_temp) %*% t(y_temp) / (tt - length(pos))
          } 
        } else {
          diag(s_endo) <- apply(matrix(y, k), 1, stats::var)
        }
      }
      if (sigma == "VAR") { # VAR model
        s_endo <- ols_sigma
      }
      s_endo <- sqrt(diag(s_endo)) # Residual standard deviations (OLS)
      
      # Endogenous variables
      if (p_endogen > 0) {
        for (r in 1:p_endogen) {
          for (l in 1:k) {
            for (j in 1:k) {
              if (l == j) {
                V[l, n_ect + (r - 1) * k + j] <- kappa1 / r^2
              } else {
                V[l, n_ect + (r - 1) * k + j] <- kappa1 * kappa2 / r^2 * s_endo[l]^2 / s_endo[j]^2
              }
            } 
          }
        } 
      }
      
      # Foreign variables
      if (k_exogen > 0) {
        p_exogen <- p_exogen - 1
        s_exo <- sqrt(apply(matrix(x[n_ect + k * p_endogen + 1:k_exogen, ], k_exogen), 1, stats::var))
        for (r in 1:(p_exogen + 1)) {
          for (l in 1:k) {
            for (j in 1:k_exogen) {
              # Note that in the loop r starts at 1, so that this is equivalent to l + 1
              if (is.null(kappa3)) {
                V[l, n_ect + k * p_endogen + (r - 1) * k_exogen + j] <- kappa1 * kappa4 * s_endo[l]^2
              } else {
                V[l, n_ect + k * p_endogen + (r - 1) * k_exogen + j] <- kappa1 * kappa3 / r^2 * s_endo[l]^2 / s_exo[j]^2 
              }
            }
          }
        } 
      }
      
      # Global variables
      if (m > 0) {
        s <- s - 1
        s_exo <- sqrt(apply(matrix(x[n_ect + k * p_endogen + k_exogen * p_exogen + 1:m,], m), 1, stats::var))
        for (r in 1:(s + 1)) {
          for (l in 1:k) {
            for (j in 1:m) {
              # Note that in the loop r starts at 1, so that this is equivalent to l + 1
              if (is.null(kappa3)) {
                V[l, n_ect + k * p_endogen + k_exogen * p_exogen + (r - 1) * m + j] <- kappa1 * kappa4 * s_endo[l]^2
              } else {
                V[l, n_ect + k * p_endogen + k_exogen * p_exogen + (r - 1) * m + j] <- kappa1 * kappa3 / r^2 * s_endo[l]^2 / s_exo[j]^2 
              }
            }
          }
        } 
      }
      
      # Restrict prior variances
      if (!is.null(max_var)) {
        if (any(stats::na.omit(c(V)) > max_var)) {
          V[which(V > max_var)] <- max_var
        } 
      }
      
      # Deterministic variables
      if (object[["model"]][["n"]] > 0){
        V[, -(1:(n_ect + k * p_endogen + k_exogen * p_exogen + m * s))] <- kappa1 * kappa4 * s_endo^2 
      }
      
      # Drop cointegration priors
      V <- V[, -(1:n_ect)]
      tot_par <- k * NCOL(object[["data"]][["train"]][["x"]])
      
      
      # Prior means
      mu <- matrix(rep(0, tot_par), k)
      mu <- matrix(mu)
      
      V <- matrix(V)
      
    } else {
      s_endo <- sqrt(diag(stats::var(t(y))))
    }
    
    
    # Structural parameters
    if (object[["model"]][["structural"]] & k > 1) {
      mu <- rbind(mu, matrix(0, k * (k - 1) / 2))
      
      V_struct <- matrix(NA, k, k)
      for (j in 1:(k - 1)) {
        V_struct[(j + 1):k, j] <- kappa1 * kappa2 * s_endo[(j + 1):k]^2 / s_endo[j]^2  
      }
      V_struct <- matrix(V_struct[lower.tri(V_struct)])
      V <- rbind(V, V_struct)
    }
    
    if (object[["model"]][["rank"]] > 0) {
      n_alpha <- k * object[["model"]][["rank"]]
      mu <- rbind(matrix(0, n_alpha), mu)
      V <- rbind(matrix(kappa1, n_alpha), V)
    }
    
    # Prior precision
    v_i <- diag(c(1 / V))
    
    result <- list("mu" = mu,
                   "v_i" = v_i)
    
    if (!is.null(object[["data"]][["train"]][["x"]])) {
      if (sigma == "AR") {
        result[["sigma_inv"]] <- matrix(0, k, k)
        diag(result[["sigma_inv"]]) = 1 / s_endo^2
      }
      if (sigma == "VAR") {
        result[["sigma_inv"]] = solve(ols_sigma)  
      }
    }
    
  }
  
  return(result)
}