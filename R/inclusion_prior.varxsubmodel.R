#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class \code{"varxsubmodel"}, usually, a result of a
#' call to \code{\link{create_varxsubmodel}}.
#' @param prob a numeric specifying the prior inclusion probability of all model parameters.
#' @param exclude_deterministics logical. If \code{TRUE} (default), the vector of the positions of
#' included variables does not include the positions of deterministic terms.
#' @param minnesota_like logical. If \code{TRUE}, the prior inclusion probabilities of the
#' parameters are calculated in a similar way as the Minnesota prior. See 'Details'.
#' @param kappa1 a numeric specifying the prior inclusion probability of
#' coefficients that correspond to own lags of endogenous variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa2 a numeric specifying the size of the prior inclusion probabilities
#' of endogenous variables, which do not correspond to own lags.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa3 a numeric specifying the size of the prior inclusion probabilities
#' of non-deterministic exogenous variables. Default is \code{NULL}, which indicates that the formula
#' for the calculation of the prior inclusion probabilities of deterministic terms
#' is used for all exogenous variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa4 a numeric specifying the size of the prior inclusion probabilities
#' of deterministic terms. Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' 
#' @details If \code{minnesota_like = TRUE}, prior inclusion probabilities \eqn{\underline{\pi}_1}
#' are calculated as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{r}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{r}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + r}} \tab for unmodelled exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables.
#' }
#' 
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selection algorithm.
#' 
#' @examples
#' 
#' # Load data
#' data("gvar2019")
#' global_data <- gvar2019[["global_data"]]
#' submodel_data <- gvar2019[["submodel_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#' 
#' # Create global model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2013:2016)
#' 
#' # Create sub-model
#' object <- create_varxsubmodel(object,
#'                               submodel = "AT",
#'                               endogen = c("y","Dp", "r"), p_endogen = 1,
#'                               exogen = c("y", "Dp"), p_exogen = 1,
#'                               global = "poil", s = 0,
#'                               deterministic = "const", seasonal = FALSE,
#'                               structural = FALSE, tvp = FALSE,
#'                               error = "wishart", varsel = "none",
#'                               iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # The previous function returns a model list. Extract the first model to proceed.
#' object <- object[[1]]
#' 
#' prior <- inclusion_prior(object)
#' 
#' @export
#' @method inclusion_prior varxsubmodel
inclusion_prior.varxsubmodel <- function(object,
                                         prob = .5,
                                         exclude_deterministics = TRUE,
                                         minnesota_like = FALSE,
                                         kappa1 = 0.8,
                                         kappa2 = 0.5,
                                         kappa3 = 0.5,
                                         kappa4 = 0.8) {
  
  if (!minnesota_like) {
    if (prob > 1 | prob < 0) {
      stop("Argument 'prob' must be between 0 and 1.")
    } 
  }
  if (minnesota_like) {
    
    if (kappa1 < 0) {
      stop("Argument 'kappa1' must not be negative.")
    }
    if (kappa1 > 1) {
      stop("Argument 'kappa1' must not be larger than 1.")
    }
    
    if (kappa2 <= 0) {
      stop("Argument 'kappa2' must not be negative.")
    }
    if (kappa2 > 1) {
      stop("Argument 'kappa2' must not be larger than 1.")
    }
    
    if (!is.null(kappa3)) {
      if (kappa3 <= 0) {
        stop("Argument 'kappa3' must not be negative.")
      } 
      if (kappa3 > 1) {
        stop("Argument 'kappa3' must not be larger than 1.")
      }
    }
    
    if (kappa4 <= 0) {
      stop("Argument 'kappa4' must not be negative.")
    }
    if (kappa4 > 1) {
      stop("Argument 'kappa4' must not be larger than 1.")
    }
  }
  
  result <- NULL
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    
    z <- object[["data"]][["train"]][["z"]]
    k <- object[["model"]][["k"]]
    tt <- nrow(object[["data"]][["train"]][["y"]])
    k_endogen <- object[["model"]][["k_endogen"]]
    p_endogen <- object[["model"]][["p_endogen"]]
    n_endogen <- k_endogen * p_endogen
    k_exogen <- object[["model"]][["k_exogen"]]
    p_exogen <- object[["model"]][["p_exogen"]]
    n_exogen <- k_exogen * (p_exogen + 1)
    m <- object[["model"]][["m_global"]]
    s <- object[["model"]][["s_global"]]
    n_global <- m * (s + 1)
    n_c <- object[["model"]][["n"]]
    
    inprior <- rep(prob, ncol(z))
    include <- 1:ncol(z)
    
    if (minnesota_like) {
      
      incl_matrix <- matrix(NA, k_endogen, n_endogen + n_exogen + n_global + n_c)
      
      if (p_endogen > 0) {
        for (i in 1:p_endogen) {
          incl_matrix[, (i - 1) * k_endogen + 1:k_endogen] <- kappa2 / i
          if (k_endogen > 1) {
            diag(incl_matrix[, (i - 1) * k_endogen + 1:k_endogen]) <- kappa1 / i 
          } else {
            incl_matrix[, (i - 1) * k_endogen + 1] <- kappa1 / i
          }
        }
      }
      
      if (k_exogen > 0) {
        incl_matrix[, n_endogen + 1:k_exogen] <- kappa3
        if (p_exogen > 0) {
          for (i in 1:p_exogen) {
            incl_matrix[, n_endogen + k_exogen + (i - 1) * k_exogen + 1:k_exogen] <- kappa3 / (1 + i)
          }
        }
      }
      
      if (m > 0) {
        incl_matrix[, n_endogen + n_exogen + 1:m] <- kappa3
        if (s > 0) {
          for (i in 1:s) {
            incl_matrix[, n_endogen + n_exogen + m + (i - 1) * m + 1:m] <- kappa3 / (1 + i)
          }
        }
      }
      
      if (n_c > 0) {
        incl_matrix[, n_endogen + n_exogen + n_global + 1:n_c] <- kappa4
      }
      
      inprior[1:(k_endogen * (n_endogen + n_exogen + n_global + n_c))] <- c(incl_matrix)
    }
    
    # Exclude deterministics from variables selection algorithm
    if (n_c > 0 & exclude_deterministics) {
      pos_det <- k_endogen * (n_endogen + n_exogen + n_global) + 1:(k_endogen * n_c)
      include <- include[-pos_det]
    }
    
    if (length(include) > 0) {
      result <- list("prior" = matrix(inprior),
                     "include" = matrix(include))
    }
  }
  
  return(result)
}