#' Summarising Country-Specific VARX Submodels of a GVAR Model
#'
#' summary method for class \code{"ctryvarest"}.
#'
#' @param object an object of class \code{"ctryvarest"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param period integer. Index of the period, for which the summary statistics should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param x an object of class \code{"summary.ctryvarest"}, usually, a result of a call to
#' \code{\link{summary.ctryvarest}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.ctryvarest} returns a list of class \code{"summary.ctryvarest"},
#' which contains the following components:
#' \item{coefficients}{A list of various summary statistics of the posterior
#' draws of the VAR coefficients.}
#' \item{sigma}{A list of various summary statistics of the posterior
#' draws of the variance-covariance matrix.}
#' \item{model}{a list containing information on the model specifications.}
#'
#' @export
summary.ctryvarest <- function(object, ci = .95, period = NULL, ...){
  
  # Get model specs ----
  names_domestic <- object[["model"]][["domestic"]][["variables"]]
  k_domestic <- length(object[["model"]][["domestic"]][["variables"]])
  p_domestic <- object[["model"]][["domestic"]][["lags"]]
  names_foreign <- object[["model"]][["foreign"]][["variables"]]
  k_foreign <- length(object[["model"]][["foreign"]][["variables"]])
  p_foreign <- object[["model"]][["foreign"]][["lags"]]
  names_global <- object[["model"]][["global"]][["variables"]]
  global <- !is.null(object[["model"]][["global"]][["variables"]])
  if (global) {
    names_global <- object[["model"]][["global"]][["variables"]]
    k_global <- length(object[["model"]][["global"]][["variables"]])
    s_global <- object[["model"]][["global"]][["lags"]] 
  } else {
    k_global <- 0
    s_global <- 0
  }
  names_deterministic <- object[["model"]][["deterministic"]]
  if (is.null(names_deterministic)) {
    k_deterministic <- 0
  } else {
    k_deterministic <- length(object[["model"]][["deterministic"]]) 
  }
  
  tt <- NROW(object[["data"]][["Y"]])
  tvp <- object[["model"]][["tvp"]]
  sv <- object[["model"]][["sv"]]
  if (tvp | sv) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
    period_long <- stats::time(object[["data"]][["Y"]])[period]
  } else {
    period_long <- NULL
  }
  
  # Obtain variable names
  x_names <- .get_regressor_names_var(object)
  dim_names <- list(names_domestic, x_names)
  
  # Non-error coefficients
  means <- NULL
  median <- NULL
  sds <- NULL
  naive_sd <- NULL
  ts_sd <- NULL
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  q_low <- NULL
  q_high <- NULL
  
  use_incl <- FALSE
  if (any(grepl("lambda", names(object[["posteriors"]])))) {
    use_incl <- TRUE
    incl <- NULL
  }
  
  vars <- c("domestic", "foreign", "global", "deterministic", "a0")
  for (i in vars) {
    if (!is.null(object[["posteriors"]][[i]])) {
      if (tvp) {
       temp <- summary(object[["posteriors"]][[i]][[period]], quantiles = c(ci_low, .5, ci_high))
      } else {
        temp <- summary(object[["posteriors"]][[i]], quantiles = c(ci_low, .5, ci_high)) 
      }
      if ("numeric" %in% class(temp[["statistics"]])) {
        means <- cbind(means, matrix(temp[["statistics"]]["Mean"], k_domestic))
        sds <- cbind(sds, matrix(temp[["statistics"]]["SD"], k_domestic))
        naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]]["Naive SE"], k_domestic))
        ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]]["Time-series SE"], k_domestic))
        q_low <- cbind(q_low, matrix(temp[["quantiles"]][1], k_domestic))
        median <- cbind(median, matrix(temp[["quantiles"]][2], k_domestic))
        q_high <- cbind(q_high, matrix(temp[["quantiles"]][3], k_domestic)) 
      } else {
        means <- cbind(means, matrix(temp[["statistics"]][, "Mean"], k_domestic))
        sds <- cbind(sds, matrix(temp[["statistics"]][, "SD"], k_domestic))
        naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]][, "Naive SE"], k_domestic))
        ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]][, "Time-series SE"], k_domestic))
        q_low <- cbind(q_low, matrix(temp[["quantiles"]][, 1], k_domestic))
        median <- cbind(median, matrix(temp[["quantiles"]][, 2], k_domestic))
        q_high <- cbind(q_high, matrix(temp[["quantiles"]][, 3], k_domestic)) 
      }
      if (use_incl) {
        var_temp <- paste0(i, "_lambda")
        if (var_temp %in% names(object[["posteriors"]])) {
          incl <- cbind(incl, matrix(colMeans(object[["posteriors"]][[var_temp]]), k_domestic))
        } else {
          incl <- cbind(incl, matrix(rep(NA_real_, ncol(object[["posteriors"]][[i]])), k_domestic))
        }
      }
    }
  }
  
  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
    if (use_incl) {
      dimnames(incl) <- dim_names
    }
  }
  
  result <- list(coefficients = list(means = means,
                                     median = median,
                                     sd = sds,
                                     naivesd = naive_sd,
                                     tssd = ts_sd,
                                     q_lower = q_low,
                                     q_upper = q_high))
  
  if (use_incl) {
    result[["coefficients"]][["lambda"]] = incl
  }
  
  dim_names <- list(dim_names[[1]], dim_names[[1]])
  
  # Error coefficients
  if (!is.null(object[["posteriors"]][["sigma"]])) {
    if (sv) {
      temp <- summary(object[["posteriors"]][["sigma"]][[period]], quantiles = c(ci_low, .5, ci_high))
    } else {
      temp <- summary(object[["posteriors"]][["sigma"]], quantiles = c(ci_low, .5, ci_high)) 
    }
    if (k_domestic == 1) {
      means <- matrix(temp[["statistics"]]["Mean"], k_domestic)
      sds <- matrix(temp[["statistics"]]["SD"], k_domestic)
      naive_sd <- matrix(temp[["statistics"]]["Naive SE"], k_domestic)
      ts_sd <- matrix(temp[["statistics"]]["Time-series SE"], k_domestic)
      q_low <- matrix(temp[["quantiles"]][1], k_domestic)
      median <- matrix(temp[["quantiles"]][2], k_domestic)
      q_high <- matrix(temp[["quantiles"]][3], k_domestic)
    } else {
      means <- matrix(temp[["statistics"]][, "Mean"], k_domestic)
      sds <- matrix(temp[["statistics"]][, "SD"], k_domestic)
      naive_sd <- matrix(temp[["statistics"]][, "Naive SE"], k_domestic)
      ts_sd <- matrix(temp[["statistics"]][, "Time-series SE"], k_domestic)
      q_low <- matrix(temp[["quantiles"]][, 1], k_domestic)
      median <- matrix(temp[["quantiles"]][, 2], k_domestic)
      q_high <- matrix(temp[["quantiles"]][, 3], k_domestic)
    }

    
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
    
    result[["sigma"]] <- list(means = means,
                              median = median,
                              sd = sds,
                              naivesd = naive_sd,
                              tssd = ts_sd,
                              q_lower = q_low,
                              q_upper = q_high)
    
    if ("sigma_lambda" %in% names(object[["posteriors"]])) {
      incl <- matrix(colMeans(object[["posteriors"]][["sigma_lambda"]]), k)
      dimnames(incl) <- dim_names
      result[["sigma"]][["lambda"]] = incl
    }
  }
  
  result[["model"]] <- object[["model"]]
  result[["model"]][["ci"]] <- paste(c(ci_low, ci_high) * 100, "%", sep = "")
  result[["model"]][["period"]] <- period_long
  
  class(result) <- c("summary.ctryvarest", "list")
  return(result)
}
