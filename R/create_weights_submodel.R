
# Weights of one sub-model for every period of its endogenous series.
#
# 'period' is either the number of periods of a rolling window or the periods
# of the weight series that constant weights are summed over. A single value is
# read as the latter when it is one of the periods of the weight series and
# cannot be the length of a window over it, as with period = 2001 for annual
# weights from 1980. It used to be read as a window of 2001 years, which failed
# with a missing value in a condition. 'name' is the name of the sub-model and
# only used in messages.
.create_weights_submodel <- function(submodel, period = NULL, name = NULL) {

  endogen <- submodel[["endogen"]]
  weights <- submodel[["weights"]]

  if (!"ts" %in% class(weights)) {
    stop("Element 'weights' must be of class 'ts'.")
  }

  label <- if (is.null(name)) "" else paste0(" of sub-model ", name)

  # Time stamps of endogenous data
  tt <- NROW(endogen)
  endogen_time <- as.numeric(stats::time(endogen))
  endogen_tsp <- stats::tsp(endogen)

  # Time stamps of weight data
  weight_time <- as.numeric(stats::time(weights))
  n_weights <- length(weight_time)

  if (length(period) == 0 || anyNA(period)) {
    stop("Argument 'period' must be given and must not contain missing values.")
  }

  use_rolling_window <- length(period) == 1 &&
    !(period > n_weights && any(abs(period - weight_time) < 1e-8))

  if (use_rolling_window) {
    if (period < 1 || period != round(period)) {
      stop("A rolling window must span a positive whole number of periods, ",
           "but argument 'period' is ", period, ".")
    }
    if (period > n_weights) {
      stop("A rolling window of ", period, " periods needs at least as many ",
           "periods of weights, but the weight series", label, " has ",
           n_weights, ", from ", weight_time[1], " to ", weight_time[n_weights],
           ". To use constant weights, give periods of the weight series.")
    }
  } else {
    missing_periods <- period[!vapply(period, function(p) any(abs(p - weight_time) < 1e-8),
                                      logical(1))]
    if (length(missing_periods) > 0) {
      stop("The weight series", label, " runs from ", weight_time[1], " to ",
           weight_time[n_weights], " and does not contain the period(s) ",
           paste0(missing_periods, collapse = ", "), " of argument 'period'.")
    }
  }

  # Determine availability of weight data per period
  if (use_rolling_window) {
    availability_matrix <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (endogen_time[i] <= weight_time[period]) {
        availability_matrix[i,] <- weight_time[1:period]
      }
      if (endogen_time[i] >= weight_time[period]) {
        if (any(floor(endogen_time[i]) == weight_time)) {
          pos_t <- which(floor(endogen_time[i]) == weight_time)
          pos_t <- (pos_t - period + 1):pos_t
        }
        # If condition is not met, the pos_t, from the last iteration will be used.
        availability_matrix[i,] <- weight_time[pos_t]
      }
    }
  } else {
    availability_matrix <- t(matrix(period, nrow = length(period), ncol = tt))
  }

  # Create a weight matrix for each period
  w <- stats::ts(matrix(NA, tt, ncol(weights)),
                 start = endogen_tsp[1], frequency = endogen_tsp[3])
  dimnames(w)[[2]] <- dimnames(weights)[[2]]

  for (i in 1:tt){
    temp <- colSums(weights[stats::time(weights) %in% as.character(availability_matrix[i,]), , drop = FALSE])
    w[i, ] <- temp / sum(temp)
  }

  return(w)

}
