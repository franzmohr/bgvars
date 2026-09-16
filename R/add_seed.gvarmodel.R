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
#' \code{\link[bvartools]{add_initial_values}} already gives every model a seed
#' drawn from R's random number generator. \code{add_seed} replaces them with
#' \code{seed}, \code{seed + 1}, \ldots, one per model, counted through the
#' sub-models in their order and through the models of each sub-model in theirs.
#' Every model then draws the same numbers whichever worker of a cluster
#' simulates it, and no two models draw the same. See
#' \code{\link[bvartools]{add_seed}} for how the seed is used.
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

# The models of all sub-models as one 'modellist', so that bvartools numbers and
# checks the seeds, and then back where they came from, with the classes of the
# sub-model lists left as they were.
.add_seed_to_submodels <- function(object, seed, ...) {
  submodels <- object[["submodels"]]
  sizes <- vapply(submodels, length, integer(1))

  flat <- unlist(lapply(submodels, unclass), recursive = FALSE, use.names = FALSE)
  if (is.null(flat)) {
    flat <- list()
  }
  class(flat) <- c("modellist", "list")
  flat <- bvartools::add_seed(flat, seed, ...)

  before <- cumsum(sizes) - sizes
  for (s in seq_along(submodels)) {
    models <- submodels[[s]]
    for (p in seq_len(sizes[s])) {
      models[[p]] <- flat[[before[s] + p]]
    }
    submodels[[s]] <- models
  }

  object[["submodels"]] <- submodels
  object
}
