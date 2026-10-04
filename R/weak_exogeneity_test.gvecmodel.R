#' @include weak_exogeneity_test.R
#' @include weak_exogeneity_test.vecxsubmodel.R
NULL

#' @rdname weak_exogeneity_test.vecxsubmodel
#' @export
#' @method weak_exogeneity_test gvecmodel
weak_exogeneity_test.gvecmodel <- function(object, ...) {
  .read_from_submodel_models(object, weak_exogeneity_test, ...)
}

#' @rdname weak_exogeneity_test.vecxsubmodel
#' @export
#' @method weak_exogeneity_test gvarmodel
weak_exogeneity_test.gvarmodel <- function(object, ...) {
  stop("The sub-models of a global VAR have no cointegrating relations, so ",
       "weak exogeneity has no content for them. The assumption the test ",
       "examines is that a country's disequilibria do not enter the equations ",
       "of its foreign variables, and a model without an error correction term ",
       "specifies no disequilibrium. Build the global model with ",
       "create_gvecmodel() to ask the question.", call. = FALSE)
}

#' @rdname weak_exogeneity_test.vecxsubmodel
#' @export
#' @method weak_exogeneity_test gvarexport
weak_exogeneity_test.gvarexport <- function(object, ..., submodels = NULL,
                                            cores = 1) {
  map_submodels(object, function(model, ...) weak_exogeneity_test(model, ...),
                ..., submodels = submodels, write = FALSE, cores = cores)
}
