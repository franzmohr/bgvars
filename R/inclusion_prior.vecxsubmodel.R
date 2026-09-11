#' Prior Inclusion Probabilities
#'
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#'
#' @param object an object of class \code{"vecxsubmodel"}, usually, a result of a
#' call to \code{\link{create_vecxsubmodel}}.
#' @param prob a numeric specifying the prior inclusion probability of all model parameters.
#' @param exclude_deterministics logical. If \code{TRUE} (default), the vector of the positions of
#' included variables does not include the positions of unrestricted deterministic terms.
#' @param minnesota_like logical. If \code{TRUE}, the prior inclusion probabilities of the
#' parameters are calculated in a similar way as the Minnesota prior. See 'Details'.
#' @param kappa1 a numeric specifying the prior inclusion probability of
#' coefficients that correspond to own lags of differenced endogenous variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa2 a numeric specifying the size of the prior inclusion probabilities
#' of differenced endogenous variables, which do not correspond to own lags.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa3 a numeric specifying the size of the prior inclusion probabilities
#' of differenced weakly exogenous and global variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa4 a numeric specifying the size of the prior inclusion probabilities
#' of unrestricted deterministic terms. Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details If \code{minnesota_like = TRUE}, prior inclusion probabilities \eqn{\underline{\pi}_1}
#' are calculated as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of differenced endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other differenced endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for differenced weakly exogenous and global variables, \cr
#' \eqn{\kappa_{4}} \tab for unrestricted deterministic terms,
#' }
#' for lag \eqn{l}.
#'
#' The loadings \eqn{\alpha} of the error correction term are never subject to variable
#' selection. They are the first \eqn{K r} coefficients of a sub-model and are dropped from
#' the vector of included positions, because the cointegration term is drawn in a step of
#' its own and the sampler rejects an \code{include} that reaches into it. The elements of
#' the cointegration matrix \eqn{\Pi = \alpha \beta^{\prime}} are therefore always in the
#' model, whatever \code{prob} is.
#'
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selection algorithm.
#'
#' @references
#'
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#'
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#'
#' @examples
#'
#' # Load data
#' data("gvar2019")
#' submodel_data <- gvar2019[["submodel_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#'
#' # Create global model
#' object <- create_gvecmodel(submodel_data = submodel_data)
#'
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2013:2016)
#'
#' # Create sub-model
#' object <- create_vecxsubmodel(object,
#'                               submodel = "AT",
#'                               endogen = c("y", "Dp", "r"), p_endogen = 2,
#'                               exogen = c("y", "Dp"), p_exogen = 1,
#'                               r = 1, const = "unrestricted",
#'                               iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' # The previous function returns a model list. Extract the first model to proceed.
#' object <- object[[1]]
#'
#' prior <- inclusion_prior(object)
#'
#' @export
#' @method inclusion_prior vecxsubmodel
inclusion_prior.vecxsubmodel <- function(object,
                                         prob = .5,
                                         exclude_deterministics = TRUE,
                                         minnesota_like = FALSE,
                                         kappa1 = 0.8,
                                         kappa2 = 0.5,
                                         kappa3 = 0.5,
                                         kappa4 = 0.8, ...) {

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
    k_endogen <- object[["model"]][["k_endogen"]]
    r <- object[["model"]][["rank"]]
    n_alpha <- k_endogen * r

    # The specification carries the lag order of the level VAR, so the model
    # contains one lag of the differences less than that. The blocks of the
    # weakly exogenous and of the global variables already count their lags,
    # the first of which is the contemporaneous one.
    p_endogen <- object[["model"]][["p_endogen"]] - 1
    k_exogen <- object[["model"]][["k_exogen"]]
    p_exogen <- object[["model"]][["p_exogen"]]
    m <- object[["model"]][["m_global"]]
    s <- object[["model"]][["s_global"]]
    n_c_unres <- object[["model"]][["n"]]

    n_endogen <- k_endogen * p_endogen
    n_exogen <- k_exogen * p_exogen
    n_global <- m * s

    inprior <- rep(prob, ncol(z))
    exclude <- NULL
    include <- 1:ncol(z)

    # Loadings ----
    # Never selected, see 'Details'.
    if (r > 0) {
      inprior[1:n_alpha] <- NA
      exclude <- append(exclude, 1:n_alpha)
    }

    # Remaining coefficients ----
    if (minnesota_like) {

      incl_matrix <- matrix(NA, k_endogen, n_endogen + n_exogen + n_global + n_c_unres)

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

      # Block i holds lag i - 1, so the contemporaneous block gets kappa3 and
      # the lag l block kappa3 / (1 + l), as for a VARX sub-model.
      if (k_exogen > 0 & p_exogen > 0) {
        for (i in 1:p_exogen) {
          incl_matrix[, n_endogen + (i - 1) * k_exogen + 1:k_exogen] <- kappa3 / i
        }
      }

      if (m > 0 & s > 0) {
        for (i in 1:s) {
          incl_matrix[, n_endogen + n_exogen + (i - 1) * m + 1:m] <- kappa3 / i
        }
      }

      if (n_c_unres > 0) {
        incl_matrix[, n_endogen + n_exogen + n_global + 1:n_c_unres] <- kappa4
      }

      inprior[n_alpha + 1:(k_endogen * (n_endogen + n_exogen + n_global + n_c_unres))] <- c(incl_matrix)
    }

    # Exclude deterministics from the variable selection algorithm
    if (n_c_unres > 0 & exclude_deterministics) {
      exclude <- append(exclude,
                        n_alpha + k_endogen * (n_endogen + n_exogen + n_global) + 1:(k_endogen * n_c_unres))
    }

    if (length(exclude) > 0) {
      include <- include[-exclude]
    }

    if (length(include) > 0) {
      result <- list("prior" = matrix(inprior),
                     "include" = matrix(include))
    }
  }

  return(result)
}
