#' @include weak_exogeneity_test.R
#' @include submodel_models.R
NULL

#' Test the Weak Exogeneity of the Foreign Variables of a Sub-Model
#'
#' Tests whether the cointegrating relations of a VEC sub-model enter the
#' equations of the variables that sub-model treats as weakly exogenous.
#'
#' @param object an object of class 'vecxsubmodel' carrying draws of the
#' cointegration space, a 'modellist' of them, a 'gvecmodel', or a 'gvarexport'.
#' @param lags the number of lags of the differenced foreign variables in the
#' marginal model. Defaults to one.
#' @param level the level the share of rejecting draws is counted at. Defaults
#' to 0.05.
#' @param ... further arguments passed to or from other methods.
#' @param submodels character vector of the sub-models the test is computed for,
#' or \code{NULL}, the default, for all of them. Only the method for a
#' 'gvarexport' takes it, since only there is a model read from a file.
#' @param cores the number of worker processes the models are read on, for the
#' method of a 'gvarexport'. Defaults to one.
#'
#' @details
#' \strong{Weak exogeneity is the assumption a global model is solved on.} A
#' sub-model is estimated conditional on its foreign variables, and
#' \code{\link{submodels_to_gvar}} may stack the conditional models into one
#' system only if those variables are long-run forcing for the country: the
#' country's disequilibria must not feed back into them. The assumption is
#' therefore a restriction on the \emph{marginal} model of the foreign
#' variables, which is never estimated, and it has to be tested separately.
#'
#' The test is the one of Harbo, Johansen, Nielsen and Rahbek (1998) as
#' Dees et al. (2007) apply it. For each weakly exogenous variable the marginal
#' model
#'
#' \deqn{\Delta x^{*}_{\ell t} = \mu_{\ell} + \sum_{j=1}^{r} \gamma_{\ell j}
#' ECM^{j}_{t-1} + \Phi_{\ell}' \Delta z_{t} + \eta_{\ell t}}
#'
#' is estimated, where \eqn{ECM^{j}_{t-1} = \beta_{j}' w_{t-1}} is the j-th
#' cointegrating relation of the sub-model and \eqn{\Delta z_{t}} collects the
#' lagged differences the conditional model itself uses -- the domestic ones,
#' the global ones where the model has a global block -- together with
#' \code{lags} lags of the differenced foreign variables. The null hypothesis is
#' \eqn{\gamma_{\ell 1} = \dots = \gamma_{\ell r} = 0}, and the statistic is the
#' F test of those \eqn{r} restrictions.
#'
#' \strong{A rejection says the sub-model should not be stacked as it stands.}
#' The variable concerned is not long-run forcing for that country: its own
#' equation adjusts to the country's disequilibrium, so treating it as given is
#' wrong and the global solution rests on a restriction the data refuse. What
#' the literature does about it is to move the variable into the endogenous
#' block of the country whose model it belongs to, or to accept the rejection
#' and say so.
#'
#' \strong{The statistic is computed at every draw of the cointegration space},
#' not once at a point estimate, because \eqn{\beta} is what the test conditions
#' on and a global model's cointegrating vectors are not estimated precisely.
#' The regression itself is least squares, as in the paper. What comes back is
#' therefore the posterior median of the statistic and of its p value, together
#' with the share of draws that reject at \code{level}, which is the honest
#' summary of a test whose regressor is itself uncertain: a variable whose
#' rejection depends on which draw of \eqn{\beta} is taken has not been shown to
#' violate the assumption.
#'
#' The sample is that of the conditional model, less the rows the lags of the
#' foreign differences cost. Degrees of freedom are taken from the rank of the
#' designs rather than from their column counts, so a collinear control --
#' which a model with a restricted trend and few lags can produce -- reduces
#' them rather than breaking the test.
#'
#' @section Rank zero and time varying coefficients:
#' At rank zero there are no relations for the foreign variables to be weakly
#' exogenous with respect to, and the method says so rather than returning a
#' test of nothing.
#'
#' A sub-model with time varying parameters carries one cointegration space per
#' period, so its error correction term is a different combination of the same
#' variables in every quarter. The method refuses such a model: the marginal
#' model above is one regression with one regressor per relation, and what the
#' corresponding object is when the relation drifts is a question this test does
#' not answer. Test the specification with constant coefficients and read the
#' answer as being about the data rather than about that one specification.
#'
#' @return A data frame with one row per weakly exogenous variable and the
#' columns
#' \item{variable}{the weakly exogenous variable whose marginal model was
#' estimated.}
#' \item{rank}{the number of cointegrating relations tested, the numerator
#' restrictions.}
#' \item{lags}{the lags of the differenced foreign variables in the marginal
#' model.}
#' \item{observations}{the length of the test sample.}
#' \item{df1, df2}{the degrees of freedom of the F test.}
#' \item{statistic}{the posterior median of the F statistic.}
#' \item{p_value}{the posterior median of its p value.}
#' \item{rejected}{the share of draws whose p value falls below \code{level}.}
#'
#' For a 'modellist' a list of such data frames, one per model. For a
#' 'gvecmodel' one list per sub-model, each holding the result of every candidate
#' model of that sub-model. For a 'gvarexport' one flat named list, named by the
#' sub-model and the file that holds the model, as in \code{AT/001}.
#'
#' @references
#'
#' Dees, S., di Mauro, F., Pesaran, M. H., & Smith, L. V. (2007). Exploring the
#' international linkages of the euro area: A global VAR analysis.
#' \emph{Journal of Applied Econometrics, 22}(1), 1--38. \doi{10.1002/jae.932}
#'
#' Harbo, I., Johansen, S., Nielsen, B., & Rahbek, A. (1998). Asymptotic
#' inference on cointegrating rank in partial systems. \emph{Journal of Business
#' & Economic Statistics, 16}(4), 388--399. \doi{10.1080/07350015.1998.10524779}
#'
#' @examples
#'
#' data("gvar2023")
#' submodel_data <- select_list_elements(gvar2023[["submodel_data"]],
#'                                       c("US", "JP", "CA"))
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = gvar2023[["global_data"]])
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data, period = 3)
#' object <- add_submodels(object, endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 2,
#'                         const = "unrestricted", trend = "restricted", r = 1,
#'                         iterations = 50, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' object <- add_priors(object, coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' tests <- weak_exogeneity_test(object)
#' tests[["US"]][[1]]
#'
#' @export
#' @method weak_exogeneity_test vecxsubmodel
weak_exogeneity_test.vecxsubmodel <- function(object, lags = 1, level = 0.05,
                                              ...) {

  if (!is.numeric(lags) || length(lags) != 1 || is.na(lags) || lags < 0 ||
        lags != round(lags)) {
    stop("Argument 'lags' must be a single non-negative integer.", call. = FALSE)
  }
  lags <- as.integer(lags)
  if (!is.numeric(level) || length(level) != 1 || is.na(level) ||
        level <= 0 || level >= 1) {
    stop("Argument 'level' must be a single number between 0 and 1.",
         call. = FALSE)
  }

  specification <- object[["model"]]
  rank <- specification[["rank"]]
  if (is.null(rank) || rank < 1) {
    stop("The sub-model has rank zero, so it has no cointegrating relation for ",
         "its foreign variables to be weakly exogenous with respect to. The ",
         "test asks whether the relations enter their equations and there are ",
         "none to enter.", call. = FALSE)
  }

  beta <- object[["posterior"]][["beta"]][["coeffs"]]
  if (is.null(beta)) {
    stop("The sub-model carries no draws of the cointegration space. Call ",
         "add_posterior_coefficients() on the global model before the test.",
         call. = FALSE)
  }
  beta <- as.matrix(beta)

  k_beta <- specification[["k_beta"]]
  if (ncol(beta) != k_beta * rank) {
    stop("The sub-model carries one cointegration space per period, which a ",
         "time varying specification does, and the marginal model of this test ",
         "is one regression with one regressor per relation. Test the ",
         "specification with constant coefficients instead.", call. = FALSE)
  }

  exogen <- object[["data"]][["original"]][["exogen"]]
  if (is.null(exogen) || ncol(exogen) < 1) {
    stop("The sub-model has no weakly exogenous variables, so there is nothing ",
         "to test.", call. = FALSE)
  }
  w <- object[["data"]][["train"]][["w"]]
  x <- object[["data"]][["train"]][["x"]]
  y <- object[["data"]][["train"]][["y"]]

  variables <- dimnames(exogen)[[2]]
  differences <- diff(exogen)

  # The test sample. The conditional model's rows are those of 'y' and 'w'; the
  # lags of the foreign differences cost the first few of them, and they are
  # dropped rather than the series being padded.
  positions <- match(round(as.numeric(stats::time(y)), 6),
                     round(as.numeric(stats::time(differences)), 6))
  if (anyNA(positions)) {
    stop("The differences of the weakly exogenous variables do not cover the ",
         "sample of the sub-model, so the marginal model cannot be aligned ",
         "with the conditional one.", call. = FALSE)
  }
  usable <- which(positions - lags >= 1L)
  if (length(usable) < 1) {
    stop("Argument 'lags' asks for more lags of the foreign differences than ",
         "the sample has before the first observation of the sub-model.",
         call. = FALSE)
  }
  rows_model <- usable
  rows_differences <- positions[usable]
  observations <- length(usable)

  dependent <- unclass(differences)[rows_differences, , drop = FALSE]
  error_correction <- unclass(w)[rows_model, , drop = FALSE]

  # The controls: the lagged differences the conditional model itself carries,
  # less the contemporaneous and lagged foreign differences, which are the
  # dependent variables of the marginal model, and less the unrestricted
  # deterministic terms, which the intercept below replaces. A global block
  # stays, as the paper has it.
  controls <- NULL
  if (!is.null(x)) {
    names_x <- dimnames(x)[[2]]
    prefixes <- paste0("d.", variables, ".")
    foreign <- vapply(names_x, function(name) any(startsWith(name, prefixes)),
                      logical(1))
    deterministic <- names_x %in% specification[["deterministic"]]
    keep <- !(foreign | deterministic)
    if (any(keep)) {
      controls <- unclass(x)[rows_model, keep, drop = FALSE]
    }
  }
  if (lags > 0) {
    lagged <- do.call(cbind, lapply(seq_len(lags), function(l) {
      unclass(differences)[rows_differences - l, , drop = FALSE]
    }))
    colnames(lagged) <- paste0("d.", rep(variables, times = lags), ".lag",
                               rep(seq_len(lags), each = length(variables)))
    controls <- if (is.null(controls)) lagged else cbind(controls, lagged)
  }

  restricted <- cbind(intercept = rep(1, observations), controls)
  qr_restricted <- qr(restricted)
  if (observations - qr_restricted[["rank"]] - rank < 1) {
    stop("The marginal model has no degrees of freedom left: ", observations,
         " observations against ", qr_restricted[["rank"]] + rank,
         " regressors. Fewer lags, or a longer sample.", call. = FALSE)
  }
  rss_restricted <- vapply(seq_along(variables), function(l) {
    sum(qr.resid(qr_restricted, dependent[, l])^2)
  }, numeric(1))

  draws <- nrow(beta)
  statistics <- matrix(NA_real_, draws, length(variables))
  p_values <- matrix(NA_real_, draws, length(variables))
  df1 <- NA_integer_
  df2 <- NA_integer_

  for (d in seq_len(draws)) {

    terms <- error_correction %*% matrix(beta[d, ], k_beta, rank)
    unrestricted <- cbind(terms, restricted)
    qr_unrestricted <- qr(unrestricted)

    df1 <- qr_unrestricted[["rank"]] - qr_restricted[["rank"]]
    df2 <- observations - qr_unrestricted[["rank"]]
    if (df1 < 1 || df2 < 1) {
      next
    }

    for (l in seq_along(variables)) {
      rss <- sum(qr.resid(qr_unrestricted, dependent[, l])^2)
      statistics[d, l] <- ((rss_restricted[l] - rss) / df1) / (rss / df2)
      p_values[d, l] <- stats::pf(statistics[d, l], df1, df2,
                                  lower.tail = FALSE)
    }
  }

  if (all(is.na(statistics))) {
    stop("The error correction terms are collinear with the controls of the ",
         "marginal model in every draw, so the restriction cannot be tested.",
         call. = FALSE)
  }

  data.frame(variable = variables,
             rank = as.integer(rank),
             lags = lags,
             observations = observations,
             df1 = df1, df2 = df2,
             statistic = apply(statistics, 2, stats::median, na.rm = TRUE),
             p_value = apply(p_values, 2, stats::median, na.rm = TRUE),
             rejected = colMeans(p_values < level, na.rm = TRUE),
             row.names = NULL,
             stringsAsFactors = FALSE)
}

#' @rdname weak_exogeneity_test.vecxsubmodel
#' @export
#' @method weak_exogeneity_test varxsubmodel
weak_exogeneity_test.varxsubmodel <- function(object, ...) {
  stop("A VARX sub-model has no cointegrating relations, so weak exogeneity ",
       "has no content for it: the assumption the test examines is that the ",
       "country's disequilibria do not enter the equations of its foreign ",
       "variables, and a model in levels or differences specifies no ",
       "disequilibrium. Build the sub-models with create_gvecmodel() to ask ",
       "the question.", call. = FALSE)
}

#' @rdname weak_exogeneity_test.vecxsubmodel
#' @export
#' @method weak_exogeneity_test modellist
weak_exogeneity_test.modellist <- function(object, ...) {
  lapply(object, weak_exogeneity_test, ...)
}
