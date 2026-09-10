# Fixtures used across the test suite.
#
# The shipped data sets are large. Tests therefore work on small subsets of a
# few sub-models and two variables, which keeps posterior simulation in the
# order of a second while still exercising the full code paths.

# Subset an object of class 'submodeldata' to a few sub-models.
#
# Weight series are reduced to the retained sub-models as well, so that the
# result still satisfies the requirements checked by .check_submodeldata.
subset_submodeldata <- function(submodel_data, submodels) {
  result <- submodel_data[submodels]
  for (i in submodels) {
    result[[i]][["weights"]] <- result[[i]][["weights"]][, submodels]
  }
  class(result) <- c("submodeldata", "list")
  result
}

# 'submodeldata' object with three countries from the gvar2023 data set.
gvar_data <- function(submodels = c("US", "JP", "CA")) {
  utils::data("gvar2023", package = "bgvars", envir = environment())
  subset_submodeldata(gvar2023[["submodel_data"]], submodels)
}

gvar_global_data <- function() {
  utils::data("gvar2023", package = "bgvars", envir = environment())
  gvar2023[["global_data"]]
}

# 'submodeldata' object with three countries from the dees2007 data set.
gvec_data <- function(submodels = c("US", "JP", "CA")) {
  utils::data("dees2007", package = "bgvars", envir = environment())
  subset_submodeldata(dees2007[["submodel_data"]], submodels)
}

gvec_global_data <- function() {
  utils::data("dees2007", package = "bgvars", envir = environment())
  dees2007[["global_data"]]
}

# A 'gvarmodel' object with weight matrices, but without sub-models.
gvar_object <- function(submodels = c("US", "JP", "CA")) {
  submodel_data <- gvar_data(submodels)
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())
  add_weight_matrices(object = object, submodel_data = submodel_data, period = 3)
}

# A 'gvecmodel' object with weight matrices, but without sub-models.
gvec_object <- function(submodels = c("US", "JP", "CA")) {
  submodel_data <- gvec_data(submodels)
  object <- create_gvecmodel(submodel_data = submodel_data,
                             global_data = gvec_global_data())
  add_weight_matrices(object = object, submodel_data = submodel_data,
                      period = 1999:2001)
}

# A fully estimated 'gvarmodel' object. Deliberately tiny: two endogenous
# variables, no lags of weakly exogenous variables and 50 draws.
gvar_estimated <- function(iterations = 50, burnin = 10, p_endogen = 1) {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = p_endogen,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          deterministic = "const",
                          iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

# A fully estimated 'gvecmodel' object.
#
# Note that no global variables are used. See test-known-issues.R.
gvec_estimated <- function(iterations = 50, burnin = 10, r = 1) {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = r,
                          iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

# A fully estimated 'gvarmodel' with control over the sub-model specification.
#
# Everything passed through '...' reaches add_submodels(), which makes it
# possible to vary the order and the selection of the endogenous variables, the
# lag orders and the use of global variables from a test.
gvar_estimated_spec <- function(..., iterations = 30, burnin = 10,
                                submodels = c("US", "JP", "CA"),
                                loglik = FALSE) {
  object <- gvar_object(submodels)
  object <- add_submodels(object, ..., iterations = iterations, burnin = burnin)
  object <- align_model_obs(object)
  object <- add_priors(object,
                       coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)
  if (loglik) {
    object <- add_posterior_loglik(object)
  }
  object
}

# An estimated 'gvarmodel' holding several candidate models per sub-model,
# together with the draws of the log-likelihood the selection criteria need.
gvar_estimated_grid <- function(iterations = 30, burnin = 10) {
  gvar_estimated_spec(endogen = c("y", "Dp"), p_endogen = 1:2,
                      exogen = c("y", "Dp"), p_exogen = 0,
                      global = "poil", s = 0,
                      deterministic = "const",
                      iterations = iterations, burnin = burnin,
                      loglik = TRUE)
}

# The Pesaran-Shin generalised impulse response of a reduced form 'bvarmodel',
# calculated from its posterior draws without using bvartools.
#
# Used to check that what submodels_to_gvar() returns really is the reduced
# form of the global model, rather than only having the right shape.
reference_girf <- function(x, impulse, response, n_ahead) {

  k <- x[["model"]][["k"]]
  p <- x[["model"]][["p"]]
  coeffs <- x[["posterior"]][["a"]][["coeffs"]]
  precision <- x[["posterior"]][["u_sigma_inv"]][["coeffs"]]

  j <- which(x[["model"]][["endogen"]] == impulse)
  i <- which(x[["model"]][["endogen"]] == response)

  one_draw <- function(draw) {
    a <- matrix(coeffs[draw, 1:(k * k * p)], k)
    sigma <- solve(matrix(precision[draw, ], k))

    phi <- list(diag(1, k))
    for (h in 1:n_ahead) {
      temp <- matrix(0, k, k)
      for (l in 1:min(h, p)) {
        temp <- temp + phi[[h - l + 1]] %*% a[, (l - 1) * k + 1:k]
      }
      phi[[h + 1]] <- temp
    }

    # One standard deviation shock to the impulse variable.
    vapply(phi, function(z) (z %*% sigma[, j])[i] / sqrt(sigma[j, j]),
           numeric(1))
  }

  vapply(seq_len(nrow(coeffs)), one_draw, numeric(n_ahead + 1))
}
