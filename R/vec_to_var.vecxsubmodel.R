#' Transform a VECX Sub-Model to a VARX Sub-Model
#'
#' Rewrites an estimated sub-model in error correction form as a sub-model in
#' levels.
#'
#' @param object an object of class 'vecxsubmodel'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The transformation is carried out here rather than delegated to
#' \code{\link[bvartools]{vec_to_var}}, because a sub-model of a global model
#' can carry a block of global variables on top of its weakly exogenous ones,
#' and the two blocks may enter with different lag orders. A 'bvecmodel'
#' describes a single block of unmodelled variables with a single lag order, so
#' it cannot express that.
#'
#' Write the sub-model as
#' \deqn{\Delta y_t = \Pi^{+} w_{t-1} + \sum_{l=1}^{p-1} \Gamma_l \Delta y_{t-l}
#' + \sum_{l=0}^{q-1} \Upsilon_l \Delta x_{t-l} + C^{ur} d^{ur}_t + u_t,}
#' where \eqn{w_{t-1}} collects the first lags of the endogenous variables, of
#' the weakly exogenous variables, of the global variables and the restricted
#' deterministic terms, and \eqn{\Pi^{+} = \alpha \beta^{+\prime}} splits into
#' \eqn{(\Pi^{d}, \Pi^{*}, \Pi^{g}, \Pi^{R})} along those blocks. Undoing the
#' differences gives the coefficients of the level representation,
#' \deqn{A_1 = I + \Pi^{d} + \Gamma_1, \quad
#' A_l = \Gamma_l - \Gamma_{l-1}, \quad A_p = -\Gamma_{p-1}}
#' for the endogenous variables, and, for each block of unmodelled variables
#' with \eqn{q} blocks of differences and level coefficient \eqn{\Pi^{x}},
#' \deqn{B_0 = \Upsilon_0, \quad B_1 = \Upsilon_1 - \Upsilon_0 + \Pi^{x}, \quad
#' B_l = \Upsilon_l - \Upsilon_{l-1}, \quad B_q = -\Upsilon_{q-1}.}
#' A block of \eqn{q} differences therefore spans the lags \eqn{0} to \eqn{q} in
#' levels. \code{p_exogen} and \code{s_global} of a VECX sub-model are the lag
#' orders of its model in levels, counted from the contemporaneous term onwards,
#' so the sub-model holds \eqn{q} = \code{p_exogen - 1} differences, and the
#' result carries the same \code{p_exogen} and \code{s_global} as the sub-model
#' it is given. A block without differences, \eqn{q = 0}, still enters through
#' the error correction term at lag one, so in levels it spans the lags \eqn{0}
#' and \eqn{1}, with \eqn{B_0 = 0} and \eqn{B_1 = \Pi^{x}}, and the result
#' counts two blocks where the sub-model counts one.
#'
#' The coefficients of the result are ordered the way an estimated VARX
#' sub-model orders them: the lags of the endogenous variables, then the weakly
#' exogenous variables from lag zero, then the global variables from lag zero,
#' then the deterministic terms. \code{\link{submodels_to_gvar}} reads exactly
#' those positions when it stacks the sub-models. A trend that entered the error
#' correction term restricted becomes an ordinary deterministic term of the
#' level representation. At cointegration rank zero there is no error correction
#' term, so restricted terms have no effect on the model and are dropped, as
#' \code{\link[bvartools]{vec_to_var}} does. \code{\link{submodels_to_gvar}}
#' gives such a sub-model coefficients of zero on the terms it lacks.
#'
#' @return An object of class 'varxsubmodel'.
#'
#' @references
#'
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series
#' analysis} (2nd ed.). Berlin: Springer.
#'
#' @examples
#'
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("US", "JP", "CA"))
#'
#' # Set up the model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 2,
#'                         const = "unrestricted", trend = "restricted", r = 1,
#'                         iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' # Rewrite one sub-model in levels
#' model <- vec_to_var(object[["submodels"]][["US"]][[1]])
#'
#' @export
#' @method vec_to_var vecxsubmodel
vec_to_var.vecxsubmodel <- function(object, ...) {

  specification <- object[["model"]]

  if (isTRUE(specification[["structural"]])) {
    stop("Structural GVEC sub-models cannot be rewritten in levels yet.\n",
         "Feel free to send a feature request.")
  }

  k <- specification[["k_endogen"]]
  p <- specification[["p_endogen"]]
  m_exogen <- specification[["k_exogen"]]
  q_exogen <- specification[["p_exogen"]]
  m_global <- specification[["m_global"]]
  q_global <- specification[["s_global"]]
  if (is.null(m_global)) {m_global <- 0L}
  if (is.null(q_global)) {q_global <- 0L}
  # The specification carries the lag orders of the model in levels, counted
  # from the contemporaneous term onwards, so a block of unmodelled variables
  # holds one difference less.
  q_exogen <- max(q_exogen - 1L, 0L)
  q_global <- max(q_global - 1L, 0L)
  n_det_ur <- specification[["n"]]
  n_det_r <- specification[["n_restricted"]]
  if (is.null(n_det_ur)) {n_det_ur <- 0L}
  if (is.null(n_det_r)) {n_det_r <- 0L}
  rank <- specification[["rank"]]
  if (is.null(rank)) {rank <- 0L}

  # The rank enters only through Pi. Without a cointegration term the levels in
  # the error correction term carry no coefficient at all.
  if (rank == 0) {n_det_r <- 0L}

  y_diff <- object[["data"]][["train"]][["y"]]
  w <- object[["data"]][["train"]][["w"]]
  x_diff <- object[["data"]][["train"]][["x"]]

  if (is.null(y_diff) || is.null(w)) {
    stop("Argument 'object' must contain the data of the VEC model in ",
         "elements 'y' and 'w'.")
  }
  if (!is.null(attr(w, "scale"))) {
    stop("The series in the error correction term are still on the scale of ",
         "'scale_error_correction'. Use 'rescale_error_correction' to put them ",
         "back on the scale of the input data before the model is transformed.")
  }

  ts_info <- stats::tsp(y_diff)
  y_diff <- as.matrix(y_diff)
  w <- as.matrix(w)
  if (is.null(x_diff)) {
    x_diff <- matrix(NA_real_, nrow(y_diff), 0)
  } else {
    x_diff <- as.matrix(x_diff)
  }

  n_gamma <- max(p - 1, 0)
  if (NCOL(w) != k + m_exogen + m_global + specification[["n_restricted"]]) {
    stop("The error correction term does not have the number of columns ",
         "implied by the model specification.")
  }
  if (NCOL(x_diff) != k * n_gamma + m_exogen * q_exogen +
                      m_global * q_global + n_det_ur) {
    stop("The non-cointegration regressors do not have the number of columns ",
         "implied by the model specification.")
  }

  # Data of the level representation ----
  #
  # The levels of the endogenous variables, recovered from the differences and
  # the levels in the error correction term, and then their lags by undoing one
  # difference at a time.
  lags_endogen <- vector("list", p)
  lags_endogen[[1]] <- w[, 1:k, drop = FALSE]
  if (p > 1) {
    for (i in 2:p) {
      lags_endogen[[i]] <- lags_endogen[[i - 1]] -
        x_diff[, (i - 2) * k + 1:k, drop = FALSE]
    }
  }

  # The same for a block of unmodelled variables. Its first lag sits in the
  # error correction term, the current value is that lag plus the current
  # difference, and each further lag undoes one more difference. A block of
  # 'n_blocks' differences therefore spans the lags 0 to 'n_blocks'.
  unmodelled_levels <- function(level_lag_1, diffs, n_vars, n_blocks) {
    if (n_vars == 0) {return(list())}
    if (n_blocks == 0) {
      # No differences at all: the variables enter at their first lag only, so
      # the current value is not identified and its coefficient is zero.
      return(list(matrix(0, nrow(level_lag_1), n_vars), level_lag_1))
    }
    result <- vector("list", n_blocks + 1)
    result[[2]] <- level_lag_1
    result[[1]] <- result[[2]] + diffs[, 1:n_vars, drop = FALSE]
    if (n_blocks > 1) {
      for (j in 2:n_blocks) {
        result[[j + 1]] <- result[[j]] -
          diffs[, (j - 1) * n_vars + 1:n_vars, drop = FALSE]
      }
    }
    result
  }

  pos <- k * n_gamma
  lags_exogen <- unmodelled_levels(
    w[, k + seq_len(m_exogen), drop = FALSE],
    x_diff[, pos + seq_len(m_exogen * q_exogen), drop = FALSE],
    m_exogen, q_exogen)

  pos <- pos + m_exogen * q_exogen
  lags_global <- unmodelled_levels(
    w[, k + m_exogen + seq_len(m_global), drop = FALSE],
    x_diff[, pos + seq_len(m_global * q_global), drop = FALSE],
    m_global, q_global)

  pos <- pos + m_global * q_global
  det_data <- NULL
  det_names <- NULL
  if (n_det_ur > 0) {
    det_data <- x_diff[, pos + 1:n_det_ur, drop = FALSE]
    det_names <- specification[["deterministic"]]
  }
  if (n_det_r > 0) {
    det_data <- cbind(det_data,
                      w[, k + m_exogen + m_global + 1:n_det_r, drop = FALSE])
    det_names <- c(det_names, specification[["deterministic_restricted"]])
  }
  if (!is.null(det_data)) {dimnames(det_data) <- list(NULL, det_names)}

  # The highest lag of a block of unmodelled variables in levels. A block of q
  # differences reaches back to lag q, and a block without differences still
  # enters at lag one through the error correction term. Counting q + 1 blocks
  # for the latter gave one name and one coefficient block too few, and the
  # rewrite failed on its own column names.
  top_exogen <- max(q_exogen, 1L)
  top_global <- max(q_global, 1L)

  # Names follow the convention of an estimated VARX sub-model.
  lag_label <- function(lags, max_lag) {
    width <- max(2, nchar(max(max_lag, 1)))
    formatC(lags, width = width, flag = "0")
  }

  x <- do.call(cbind, c(lags_endogen, lags_exogen, lags_global))
  x_names <- c(paste0(rep(specification[["endogen"]], times = p), ".",
                      rep(lag_label(1:p, p), each = k)))
  if (m_exogen > 0) {
    x_names <- c(x_names,
                 paste0(rep(specification[["exogen"]], times = top_exogen + 1),
                        ".l", rep(lag_label(0:top_exogen, top_exogen + 1),
                                  each = m_exogen)))
  }
  if (m_global > 0) {
    x_names <- c(x_names,
                 paste0(rep(specification[["global"]], times = top_global + 1),
                        ".l", rep(lag_label(0:top_global, top_global + 1),
                                  each = m_global)))
  }
  if (!is.null(det_data)) {
    x <- cbind(x, det_data)
    x_names <- c(x_names, det_names)
  }
  dimnames(x) <- list(NULL, x_names)

  y <- y_diff + w[, 1:k, drop = FALSE]
  dimnames(y) <- list(NULL, specification[["endogen"]])

  y <- stats::ts(y, start = ts_info[1], frequency = ts_info[3],
                 class = c("mts", "ts", "matrix"))
  x <- stats::ts(x, start = ts_info[1], frequency = ts_info[3],
                 class = c("mts", "ts", "matrix"))

  # Coefficients of the level representation ----
  posterior <- object[["posterior"]]
  if (!is.null(posterior)) {

    if (isTRUE(specification[["tvp"]])) {
      stop("Posterior draws of models with time varying parameters cannot be ",
           "transformed.")
    }

    coeffs <- posterior[["a"]][["coeffs"]]
    mcpar <- attr(posterior[["u_sigma_inv"]][["coeffs"]], "mcpar")
    draws <- nrow(coeffs)
    k_ect <- NCOL(w)
    n_x <- NCOL(x_diff)

    if (NCOL(coeffs) != k * (rank + n_x)) {
      stop("The posterior draws do not have the number of columns implied by ",
           "the model specification.")
    }

    # Blocks of the level representation, in the order an estimated VARX
    # sub-model stores them.
    n_out <- k * p + m_exogen * (top_exogen + 1) + m_global * (top_global + 1) +
      n_det_ur + n_det_r
    result_coeffs <- matrix(NA_real_, draws, k * n_out)

    beta <- posterior[["beta"]][["coeffs"]]

    # One block of unmodelled variables: turns the coefficients of its
    # differences and the columns of Pi that belong to it into the coefficients
    # of its levels at the lags 0 to 'n_blocks'.
    unmodelled_coefficients <- function(upsilon, pi_x, n_vars, n_blocks) {
      if (n_vars == 0) {return(NULL)}
      out <- matrix(0, k, n_vars * (n_blocks + 1))
      if (n_blocks == 0) {
        out[, 1:n_vars] <- 0
        return(cbind(matrix(0, k, n_vars), pi_x))
      }
      block <- function(j) {upsilon[, (j - 1) * n_vars + 1:n_vars, drop = FALSE]}
      out[, 1:n_vars] <- block(1)
      for (j in 1:n_blocks) {
        current <- if (j < n_blocks) block(j + 1) else matrix(0, k, n_vars)
        value <- current - block(j)
        if (j == 1) {value <- value + pi_x}
        out[, j * n_vars + 1:n_vars] <- value
      }
      out
    }

    for (draw in seq_len(draws)) {

      a <- matrix(coeffs[draw, ], k)
      pi_plus <- matrix(0, k, k_ect)
      if (rank > 0) {
        alpha <- a[, 1:rank, drop = FALSE]
        b <- matrix(beta[draw, ], k_ect, rank)
        pi_plus <- alpha %*% t(b)
      }
      a <- a[, rank + seq_len(n_x), drop = FALSE]

      # Endogenous variables
      out_endogen <- matrix(0, k, k * p)
      gamma <- function(j) {a[, (j - 1) * k + 1:k, drop = FALSE]}
      out_endogen[, 1:k] <- diag(1, k) + pi_plus[, 1:k, drop = FALSE] +
        (if (n_gamma > 0) gamma(1) else 0)
      if (p > 1) {
        for (j in 2:p) {
          current <- if (j <= n_gamma) gamma(j) else matrix(0, k, k)
          out_endogen[, (j - 1) * k + 1:k] <- current - gamma(j - 1)
        }
      }

      offset <- k * n_gamma
      out_exogen <- unmodelled_coefficients(
        a[, offset + seq_len(m_exogen * q_exogen), drop = FALSE],
        pi_plus[, k + seq_len(m_exogen), drop = FALSE], m_exogen, q_exogen)

      offset <- offset + m_exogen * q_exogen
      out_global <- unmodelled_coefficients(
        a[, offset + seq_len(m_global * q_global), drop = FALSE],
        pi_plus[, k + m_exogen + seq_len(m_global), drop = FALSE],
        m_global, q_global)

      offset <- offset + m_global * q_global
      out_det <- NULL
      if (n_det_ur > 0) {
        out_det <- a[, offset + 1:n_det_ur, drop = FALSE]
      }
      if (n_det_r > 0) {
        out_det <- cbind(out_det,
                         pi_plus[, k + m_exogen + m_global + 1:n_det_r,
                                 drop = FALSE])
      }

      result_coeffs[draw, ] <- cbind(out_endogen, out_exogen, out_global,
                                     out_det)
    }

    dimnames(result_coeffs) <- list(NULL, NULL)
    coeffs <- coda::mcmc(result_coeffs)
    if (!is.null(mcpar)) {
      attr(coeffs, "mcpar") <- mcpar
    }
    posterior[["a"]] <- list("coeffs" = coeffs)
    posterior[["beta"]] <- NULL
  }

  # Specification of the level representation ----
  model <- specification
  model[["type"]] <- if (k == 1) "ARX" else "VARX"
  model[["k"]] <- k
  model[["k_endogen"]] <- k
  model[["p"]] <- p
  model[["p_endogen"]] <- p
  model[["k_exogen"]] <- m_exogen
  model[["p_exogen"]] <- if (m_exogen > 0) top_exogen + 1L else 0L
  model[["m_global"]] <- m_global
  model[["s_global"]] <- if (m_global > 0) top_global + 1L else 0L
  model[["m"]] <- model[["k_exogen"]] * model[["p_exogen"]] +
    model[["m_global"]] * model[["s_global"]]
  model[["s"]] <- model[["p_exogen"]]
  model[["n"]] <- n_det_ur + n_det_r
  model[["deterministic"]] <- det_names
  model[["n_restricted"]] <- 0L
  model[["deterministic_restricted"]] <- NULL
  model[["k_beta"]] <- NULL

  # The rank is a property of the error correction form and says nothing about
  # the model in levels, but it is what a sub-model of a GVEC model is chosen
  # by, so it is kept for reference.
  model[["rank"]] <- specification[["rank"]]

  result <- object
  result[["model"]] <- model
  result[["data"]][["train"]] <- list("y" = y, "x" = x)
  result[["posterior"]] <- posterior

  class(result) <- c("varxsubmodel", "bvarmodel", "list")

  return(result)
}
