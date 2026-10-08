#' Seed of the Posterior Simulation
#'
#' Sets the seeds with which the posteriors of the models of a global model are
#' simulated.
#'
#' @param object an object of class 'gvarmodel' or 'gvecmodel' with sub-models.
#' @param seed a non-negative whole number no larger than
#' \code{.Machine$integer.max}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Calling this function is optional:
#' \code{\link[bayests]{add_initial_values}} already gives every model a seed
#' drawn from R's random number generator. \code{add_seed} replaces them with
#' \code{seed}, \code{seed + 1}, \ldots, one per model, counted through the
#' sub-models in their order and through the models of each sub-model in theirs.
#' Every model then draws the same numbers whichever worker of a cluster
#' simulates it, and no two models draw the same. See
#' \code{\link[bayests]{add_seed}} for how the seed is used.
#'
#' @return \code{object} with the seeds set.
#'
#' @export
#' @method add_seed gvarmodel
add_seed.gvarmodel <- function(object, seed, ...) {
  .add_seed_to_submodels(object, seed, ...)
}

#' @rdname add_seed.gvarmodel
#' @export
#' @method add_seed gvecmodel
add_seed.gvecmodel <- function(object, seed, ...) {
  .add_seed_to_submodels(object, seed, ...)
}

# All models at once, so that bvartools numbers and checks the seeds.
.add_seed_to_submodels <- function(object, seed, ...) {
  flat <- .submodel_models(object)
  models <- bayests::add_seed(flat[["models"]], seed, ...)
  .set_submodel_models(object, models, flat[["sizes"]])
}
