# The deterministic terms of a global model, on the time axis of its data.
#
# Every sub-model takes its deterministic terms from this one series rather
# than building its own. A trend that each sub-model counts from its own first
# observation is not the same regressor across sub-models: a sub-model with
# fewer lags keeps more of the early observations, so its trend is shifted
# against that of a sub-model with more lags. Each of them is internally
# consistent -- a shift of the trend is absorbed by the constant -- but the
# global model has a single trend regressor, and the sub-models have to agree
# on what it is before they can be stacked.
#
# Defining the terms once, on the time axis of the global data, makes them
# agree by construction, whatever sample and lag order a sub-model ends up
# with.
#
# The full menu is built here, because which of the terms a sub-model uses is
# only decided later, in add_submodels(). A sub-model picks the columns it
# wants by name.
.global_deterministic <- function(endogen) {

  tsp_endogen <- stats::tsp(endogen)
  tt <- nrow(endogen)

  result <- cbind(rep(1, tt), seq_len(tt))
  names_result <- c("const", "trend")

  # Seasonal dummies, likewise anchored to the global time axis: which quarter
  # or month an observation falls in is a property of the data, not of the
  # sample a sub-model happens to be estimated on.
  frequency <- stats::frequency(endogen)
  if (frequency > 1) {
    cycle <- stats::cycle(endogen)
    for (i in 1:(frequency - 1)) {
      result <- cbind(result, as.numeric(cycle == i))
      names_result <- c(names_result, paste0("season.", i))
    }
  }

  result <- stats::ts(result, start = tsp_endogen[1], frequency = tsp_endogen[3],
                      class = c("mts", "ts", "matrix"))
  dimnames(result) <- list(NULL, names_result)

  return(result)
}

# The deterministic terms one sub-model uses, taken out of the global series.
#
# 'names_wanted' are the terms the sub-model uses, in the order it wants them.
# The result keeps the time axis of the global model, so that cbind() lines it
# up with the regressors of the sub-model.
.submodel_deterministic <- function(object, names_wanted, start = NULL, end = NULL) {

  deterministic <- object[["global"]][["deterministic"]]

  if (is.null(deterministic)) {
    stop("Argument 'object' does not contain deterministic terms. It was ",
         "created by a version of this package that did not build them.")
  }

  absent <- setdiff(names_wanted, dimnames(deterministic)[[2]])
  if (length(absent) > 0) {
    stop("The global model does not provide the deterministic term(s) ",
         paste0(absent, collapse = ", "), ".")
  }

  # A caller that has already trimmed its sample asks for that window, so that
  # what comes back lines up row by row with what it has.
  if (!is.null(start) || !is.null(end)) {
    deterministic <- stats::window(deterministic, start = start, end = end)
  }

  tsp_deterministic <- stats::tsp(deterministic)

  result <- stats::ts(deterministic[, names_wanted, drop = FALSE],
                      start = tsp_deterministic[1],
                      frequency = tsp_deterministic[3],
                      class = c("mts", "ts", "matrix"))
  dimnames(result) <- list(NULL, names_wanted)

  return(result)
}
