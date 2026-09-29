# The models of all sub-models of a global model as one bvartools 'modellist',
# and back into their sub-model lists, whose classes and names are kept.
#
# A global model holds its models in one list per sub-model. bvartools' methods
# for a 'modellist' seed, simulate and evaluate the models of a list in order and,
# with 'cores' above 1, on one cluster for the whole list, so handing them all
# models at once numbers the seeds through the sub-models and shares the models
# out over the workers regardless of which sub-model they belong to.

.submodel_models <- function(object) {
  submodels <- object[["submodels"]]
  sizes <- vapply(submodels, length, integer(1))
  models <- unlist(lapply(submodels, unclass), recursive = FALSE, use.names = FALSE)
  if (is.null(models)) {
    models <- list()
  }
  class(models) <- c("modellist", "list")
  list(models = models, sizes = sizes)
}

.set_submodel_models <- function(object, models, sizes) {
  submodels <- object[["submodels"]]
  before <- cumsum(sizes) - sizes
  for (s in seq_along(submodels)) {
    list_of_models <- submodels[[s]]
    for (p in seq_len(sizes[s])) {
      list_of_models[[p]] <- models[[before[s] + p]]
    }
    submodels[[s]] <- list_of_models
  }
  object[["submodels"]] <- submodels
  object
}

# Applies the bvartools generic 'generic' to all models of 'object' at once.
.apply_to_submodel_models <- function(object, generic, ..., cores = 1) {
  flat <- .submodel_models(object)
  models <- generic(flat[["models"]], ..., cores = cores)
  .set_submodel_models(object, models, flat[["sizes"]])
}

# The same for a generic that reads rather than changes a model: the results are
# grouped back into one list per sub-model, keeping the names and the order of
# the sub-models, and the global model is not touched. A diagnostic over a global
# model returns its results, not the model, so it cannot use the function above.
# Unlike the function above it passes no 'cores': the reading generics map over
# a 'modellist' in order and take no cluster, so offering one here would
# advertise an argument that does nothing.
.read_from_submodel_models <- function(object, generic, ...) {
  flat <- .submodel_models(object)
  results <- generic(flat[["models"]], ...)
  sizes <- flat[["sizes"]]
  before <- cumsum(sizes) - sizes
  out <- lapply(seq_along(sizes), function(s) {
    if (sizes[s] == 0) {
      return(list())
    }
    results[before[s] + seq_len(sizes[s])]
  })
  names(out) <- names(object[["submodels"]])
  out
}
