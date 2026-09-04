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
