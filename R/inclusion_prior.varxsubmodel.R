#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a
#' call to \code{\link{create_var_model}}.
#' @param prob a numeric specifying the prior inclusion probability of all model parameters.
#' @param exclude_deterministics logical. If \code{TRUE} (default), the vector of the positions of
#' included variables does not include the positions of deterministic terms.
#' @param minnesota_like logical. If \code{TRUE}, the prior inclusion probabilities of the
#' parameters are calculated in a similar way as the Minnesota prior. See 'Details'.
#' @param kappa a numeric vector of four elements containing the prior inclusion probabilities
#' of coefficients that correspond to own lags of endogenous variables, to endogenous variables,
#' which do not correspond to own lags, to exogenous variables and deterministic terms, respectively.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' 
#' @details If \code{minnesota_like = TRUE}, prior inclusion probabilities \eqn{\underline{\pi}_1}
#' are calculated as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{r}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{r}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + r}} \tab for foreign and global variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{r} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' third and forth element in \code{kappa}, respectively.
#' 
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selection algorithm.
#' 
#' @examples
#' 
#' # Prepare data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model input
#' object <- create_var_model(e1)
#' 
#' # Obtain inclusion prior
#' incl <- inclusion_prior(object)
#' 
#' @export
inclusion_prior.varxsubmodel <- function(object, prob = .5, exclude_deterministics = TRUE,
                                      minnesota_like = FALSE, kappa = c(0.8, 0.5, 0.5, .8)) {
  
  if (!minnesota_like) {
    if (prob > 1 | prob < 0) {
      stop("Argument 'prob' must be between 0 and 1.")
    } 
  }
  if (minnesota_like) {
    if (any(kappa > 1) | any(kappa < 0)) {
      stop("Argument 'kappa' may only contain values between 0 and 1.")
    } 
  }
  
  result <- NULL
  if (!is.null(object[["data"]][["z"]])) {
    
    y <- t(object$data$y)
    z <- object[["data"]][["z"]]
    tt <- NCOL(y)
    k_domestic <- object$model$k_domestic
    p_domestic <- object$model$p_domestic
    n_domestic <- k_domestic * p_domestic
    k_foreign <- object$model$k_foreign
    p_foreign <- object$model$p_foreign
    n_foreign <- k_foreign * (p_foreign + 1)
    m <- object$model$m
    s <- object$model$s
    n_global <- m * (s + 1)
    n_c <- object$model$n
    
    inprior <- rep(prob, ncol(z))
    include <- 1:ncol(z)
    
    if (minnesota_like & !is.null(object$data$x)) {
      
      incl_matrix <- matrix(NA, k_domestic, n_domestic + n_foreign + n_global + n_c)
      
      if (p_domestic > 0) {
        for (i in 1:p_domestic) {
          incl_matrix[, (i - 1) * k_domestic + 1:k_domestic] <- kappa[2] / i
          if (k_domestic > 1) {
            diag(incl_matrix[, (i - 1) * k_domestic + 1:k_domestic]) <- kappa[1] / i 
          } else {
            incl_matrix[, (i - 1) * k_domestic + 1] <- kappa[1] / i
          }
        }
      }
      
      if (k_foreign > 0) {
        incl_matrix[, n_domestic + 1:k_foreign] <- kappa[3]
        if (p_foreign > 0) {
          for (i in 1:p_foreign) {
            incl_matrix[, n_domestic + k_foreign + (i - 1) * k_foreign + 1:k_foreign] <- kappa[3] / (1 + i)
          }
        }
      }
      
      if (m > 0) {
        incl_matrix[, n_domestic + n_foreign + 1:m] <- kappa[3]
        if (s > 0) {
          for (i in 1:s) {
            incl_matrix[, n_domestic + n_foreign + m + (i - 1) * m + 1:m] <- kappa[3] / (1 + i)
          }
        }
      }
      
      if (n_c > 0) {
        incl_matrix[, n_domestic + n_foreign + n_global + 1:n_c] <- kappa[4]
      }
      
      inprior[1:(k_domestic * (n_domestic + n_foreign + n_global + n_c))] <- c(incl_matrix)
    }
    
    # Exclude deterministics from variables selection algorithm
    if (n_c > 0 & exclude_deterministics) {
      pos_det <- k_domestic * (n_domestic + n_foreign + n_global) + 1:(k_domestic * n_c)
      include <- include[-pos_det]
    }
    
    if (length(include) > 0) {
      result <- list("prior" = matrix(inprior),
                     "include" = matrix(include))
    }
  }
  
  return(result)
}