#' Estimation Steps on an Exported Global Model
#'
#' The steps of an estimation, applied to the models of an export one at a time
#' rather than to a global model in memory. Each reads a model, changes it,
#' writes it back and drops it, so a step costs one model per worker.
#'
#' @param object an object of class 'gvarexport', from
#' \code{\link{read_export}}.
#' @param x the same, for the method of \code{thin}, whose generic calls its
#' first argument \code{x}.
#' @param seed an integer, the seed of the first model. The models are numbered
#' through in the order of the manifest, as they are for a global model in
#' memory.
#' @param cores the number of worker processes. Defaults to one.
#' @param export character vector of names of objects the workers need, for
#' instance the sampler that \code{posterior_function} refers to. See
#' \code{\link{map_submodels}}.
#' @param ... further arguments passed to the method for a single model.
#'
#' @details
#' The methods do what their counterparts for a 'gvarmodel' or a 'gvecmodel' do,
#' model by model, and they take the same arguments. Everything that concerns
#' the global model as a whole -- the global data, the index and the weight
#' matrices -- is written once by \code{\link[bayests]{write_to_hdf5}} and is
#' not touched again.
#'
#' \code{scale_error_correction} and \code{rescale_error_correction} carry the
#' starting values of \eqn{\beta} along with the transformation, as the method
#' for a 'gvecmodel' does.
#'
#' The models of a sub-model are drawn with a seed of their own, so the draws do
#' not depend on how many workers a step runs on, nor on whether the models were
#' estimated in one run or in several.
#'
#' A step takes \code{submodels}, which \code{\link{map_submodels}} understands:
#' a character vector of the sub-models it is applied to, the rest of the export
#' being left as it is. An estimation of many expensive models is run sub-model
#' by sub-model that way, and one that was interrupted is picked up where it
#' stopped.
#'
#' @return The object in \code{object}, with its manifest brought up to date.
#'
#' @examples
#'
#' data("gvar2023")
#' submodel_data <- select_list_elements(gvar2023[["submodel_data"]],
#'                                       c("AT", "DE", "US"))
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = gvar2023[["global_data"]])
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data, period = 3)
#' object <- add_submodels(object, endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 1,
#'                         deterministic = "const",
#'                         iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' folder <- file.path(tempdir(), "bgvars-example-steps")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
#'
#' export <- read_export(folder)
#' export <- add_priors(export, coef = list(v_i = 1),
#'                      sigma = list(df = 3, scale = 1e-8))
#' export <- add_initial_values(export)
#' export <- add_seed(export, 20260919)
#' export <- add_posterior_coefficients(export)
#' export <- add_posterior_loglik(export)
#'
#' @name export_steps
#' @family export
NULL

#' @rdname export_steps
#' @export
#' @method add_priors gvarexport
add_priors.gvarexport <- function(object, ..., cores = 1) {
  map_submodels(object, function(model, ...) add_priors(model, ...), ...,
                cores = cores)
}

#' @rdname export_steps
#' @export
#' @method add_initial_values gvarexport
add_initial_values.gvarexport <- function(object, ..., cores = 1) {
  map_submodels(object, function(model, ...) add_initial_values(model, ...), ...,
                cores = cores)
}

#' @rdname export_steps
#' @export
#' @method add_seed gvarexport
add_seed.gvarexport <- function(object, seed, ..., cores = 1) {
  # The seeds are worked out here rather than on the workers: what a worker
  # runs is sent to it as it stands, so it must not call anything the package
  # installed there may not have.
  seeds <- .model_seed(seed, seq_len(nrow(object[["manifest"]])) - 1L)
  map_submodels(object, function(model, index, seeds, ...) {
    bayests::add_seed(model, seeds[[index]], ...)
  }, seeds = seeds, ..., cores = cores)
}

#' @rdname export_steps
#' @export
#' @method add_posterior_coefficients gvarexport
add_posterior_coefficients.gvarexport <- function(object, ..., cores = 1,
                                                  export = NULL) {
  map_submodels(object, function(model, ...) add_posterior_coefficients(model, ...),
                ..., cores = cores, export = export)
}

#' @rdname export_steps
#' @export
#' @method add_posterior_loglik gvarexport
add_posterior_loglik.gvarexport <- function(object, ..., cores = 1) {
  map_submodels(object, function(model, ...) add_posterior_loglik(model, ...), ...,
                cores = cores)
}

#' @rdname export_steps
#' @export
#' @method thin gvarexport
thin.gvarexport <- function(x, ..., cores = 1) {
  map_submodels(x, function(model, ...) thin(model, ...), ..., cores = cores)
}

#' @rdname export_steps
#' @export
#' @method scale_error_correction gvarexport
scale_error_correction.gvarexport <- function(object, ..., cores = 1) {
  map_submodels(object, function(model, transform, ...) {
    model <- scale_error_correction(model, ...)
    factors <- attr(model[["data"]][["train"]][["w"]], "scale")
    if (is.null(factors)) {
      return(model)
    }
    transform(model, factors)
  }, transform = .transform_initial_beta, ..., cores = cores)
}

#' @rdname export_steps
#' @export
#' @method rescale_error_correction gvarexport
rescale_error_correction.gvarexport <- function(object, ..., cores = 1) {
  map_submodels(object, function(model, transform, ...) {
    factors <- attr(model[["data"]][["train"]][["w"]], "scale")
    model <- rescale_error_correction(model, ...)
    if (is.null(factors)) {
      return(model)
    }
    transform(model, 1 / factors)
  }, transform = .transform_initial_beta, ..., cores = cores)
}

# The seed of the model at position 'offset', counted from zero, vectorised over
# the offsets. The rule is the one bvartools' add_seed() uses for a list of
# models, so an export and a global model in memory hand out the same seeds and
# draw the same posteriors.
.model_seed <- function(seed, offset) {
  seed <- as.numeric(seed)
  if (length(seed) != 1 || is.na(seed) || seed %% 1 != 0) {
    stop("Argument 'seed' must be a single whole number.")
  }
  as.integer((seed + offset) %% (as.numeric(.Machine$integer.max) + 1))
}
