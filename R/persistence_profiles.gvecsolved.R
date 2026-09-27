#' Persistence Profiles of a Solved Global Model
#'
#' Calculates the persistence profiles of the cointegrating relations of a
#' global model solved from VEC sub-models.
#'
#' @param object an object of class 'gvecsolved', the result of a call to
#' \code{\link{submodels_to_gvar}} on sub-models of a rank of at least one.
#' @param n_ahead the number of horizons the profile is followed over.
#' @param ci a numeric between 0 and 1 specifying the probability of the
#' credible band. Defaults to 0.95.
#' @param keep_draws logical. If \code{FALSE} (default) the draws are summarised
#' by their median and the band; if \code{TRUE} they are returned as they are.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The persistence profile of Pesaran and Shin (1996) is the response of a
#' cointegrating relation to a system-wide shock, scaled to one on impact. It
#' falls to zero for a relation that is cointegrating, and the speed at which it
#' does so is the speed of convergence to equilibrium. See
#' \code{\link[bvartools]{persistence_profiles.bvecmodel}} for the expression.
#'
#' \strong{This is where the profile of a sub-model's relation belongs.} A
#' cointegrating vector of a VECX sub-model spans that country's own variables
#' and its foreign ones, and the foreign ones are not endogenous to it, so the
#' method for a single sub-model can follow only the part of the relation that
#' model governs and warns as much. In the solved global model every variable of
#' every relation is endogenous, and the relation can be followed whole.
#'
#' The vectors are those \code{\link{submodels_to_gvar}} carried over, each
#' mapped into the coordinates of the global model through the weight matrix of
#' its sub-model. There are \code{sum(r_i)} of them, which is also what the
#' count of common stochastic trends, \code{k - sum(r_i)}, is taken against.
#'
#' A global model solved at one period -- which is what a model with time
#' varying coefficients gives -- has the profiles of that period. They are a
#' statement about the system as it stood then, not about the sample.
#'
#' @return A list with one element per cointegrating relation, named for the
#' sub-model it came from, each a matrix with one row per horizon and, unless
#' \code{keep_draws} is \code{TRUE}, the columns \code{median}, \code{lower} and
#' \code{upper}.
#'
#' @references
#'
#' Dees, S., di Mauro, F., Pesaran, M. H., & Smith, L. V. (2007). Exploring the
#' international linkages of the euro area: A global VAR analysis.
#' \emph{Journal of Applied Econometrics, 22}(1), 1--38. \doi{10.1002/jae.932}
#'
#' Pesaran, M. H., & Shin, Y. (1996). Cointegration and speed of convergence to
#' equilibrium. \emph{Journal of Econometrics, 71}(1-2), 117--143.
#' \doi{10.1016/0304-4076(94)01697-6}
#'
#' @export
#' @method persistence_profiles gvecsolved
persistence_profiles.gvecsolved <- function(object, n_ahead = 20, ci = 0.95,
                                            keep_draws = FALSE, ...) {

  if (ci <= 0 || ci >= 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }

  beta <- object[["posterior"]][["beta"]][["coeffs"]]
  if (is.null(beta)) {
    stop("The solved model carries no cointegrating vectors. It was stacked ",
         "from sub-models of rank zero, or from VAR sub-models, and a model ",
         "without a long-run relation has no persistence profile.", call. = FALSE)
  }

  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  a <- as.matrix(object[["posterior"]][["a"]][["coeffs"]])
  precision <- as.matrix(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  beta <- as.matrix(beta)
  draws <- nrow(a)
  relations <- ncol(beta) / k

  profiles <- array(NA_real_, c(draws, n_ahead + 1, relations))

  for (d in seq_len(draws)) {

    coefficients <- matrix(a[d, ], k)
    covariance <- solve(matrix(precision[d, ], k))
    vectors <- matrix(beta[d, ], k, relations)

    psi <- vector("list", n_ahead + 1)
    psi[[1]] <- diag(k)
    for (h in seq_len(n_ahead)) {
      total <- matrix(0, k, k)
      for (j in seq_len(min(h, p))) {
        total <- total +
          coefficients[, (j - 1) * k + seq_len(k), drop = FALSE] %*% psi[[h - j + 1]]
      }
      psi[[h + 1]] <- total
    }

    for (j in seq_len(relations)) {
      b <- vectors[, j, drop = FALSE]
      scale <- as.numeric(crossprod(b, covariance %*% b))
      for (h in seq_len(n_ahead + 1)) {
        m <- crossprod(b, psi[[h]])
        profiles[d, h, j] <- as.numeric(m %*% covariance %*% t(m)) / scale
      }
    }
  }

  labels <- attr(object[["posterior"]][["beta"]][["coeffs"]], "relations")
  if (is.null(labels) || length(labels) != relations) {
    labels <- paste0("relation_", seq_len(relations))
  }

  low <- (1 - ci) / 2
  result <- lapply(seq_len(relations), function(j) {
    x <- matrix(profiles[, , j], nrow = draws)
    if (keep_draws) {
      colnames(x) <- 0:n_ahead
      return(x)
    }
    out <- cbind(median = apply(x, 2, stats::median),
                 lower = apply(x, 2, stats::quantile, probs = low),
                 upper = apply(x, 2, stats::quantile, probs = 1 - low))
    rownames(out) <- 0:n_ahead
    out
  })
  names(result) <- labels
  class(result) <- append("bvecpp", class(result))
  result
}
