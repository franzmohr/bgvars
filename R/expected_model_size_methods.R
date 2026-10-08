#' Expected Size of a Global Model
#'
#' Calculates how much memory, or disk space, a global model will take once the
#' posterior draws of all its sub-models are complete, from their specifications
#' alone, so that it can be called before a long estimation is started.
#'
#' @param object an object of class 'gvarmodel' or 'gvecmodel' with sub-models,
#' or an export of one of class 'gvarexport', from \code{\link{read_export}}.
#' @param ... further arguments passed to the method of \pkg{bvartools} for a
#' single model, such as \code{chains}.
#'
#' @details Every model of every sub-model is sized by
#' \code{\link[bayests]{expected_model_size}}, which describes how. A global
#' model holds all of them at once, so its size is their sum, and \strong{a
#' sub-model with time-varying parameters grows with the number of periods, the
#' draws and the number of models tried for it}. That is what an export is for:
#' it keeps the models in files and estimates them one at a time, so what it
#' needs is the disk space of the whole and the memory of one model per worker.
#'
#' For an export the models are read from their files without their draws, so
#' sizing one costs little even once it has been estimated.
#'
#' \code{\link[bayests]{add_posterior_coefficients}} warns before it starts,
#' for a global model or an export alike, when the size exceeds
#' \code{options(bayests.size_warning)}.
#'
#' @return A data frame of class 'modelsize', as described on the page of
#' \code{\link[bvartools]{expected_model_size.bvarmodel}}, with one block of rows
#' per model. The column \code{model} names each model by its sub-model and its
#' position among the models of that sub-model, such as \code{"US/2"}. The first
#' row, \code{"global"}, is what the object holds beside its sub-models: the
#' global data and the weight matrices, or for an export the file
#' \code{model.h5} that holds them.
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
#' object <- add_submodels(object, endogen = c("y", "Dp"), p_endogen = 1:2,
#'                         exogen = c("y", "Dp"), p_exogen = 1,
#'                         deterministic = "const", tvp = TRUE,
#'                         iterations = 5000, burnin = 1000)
#' expected_model_size(object)
#'
#' @export
#' @method expected_model_size gvarmodel
expected_model_size.gvarmodel <- function(object, ...) {
  .global_model_size(object, ...)
}

#' @rdname expected_model_size.gvarmodel
#' @export
#' @method expected_model_size gvecmodel
expected_model_size.gvecmodel <- function(object, ...) {
  .global_model_size(object, ...)
}

#' @rdname expected_model_size.gvarmodel
#' @export
#' @method expected_model_size gvarexport
expected_model_size.gvarexport <- function(object, ...) {

  manifest <- object[["manifest"]]
  models <- lapply(file.path(object[["folder"]], manifest[["file"]]), function(file) {
    bayests::read_model_from_hdf5(filename = file, draws = integer(0))
  })
  names(models) <- .submodel_model_labels(manifest[["submodel"]])
  class(models) <- c("modellist", "list")

  global <- file.size(file.path(object[["folder"]], "model.h5"))
  .with_global_row(bayests::expected_model_size(models, ...), global,
                   "model.h5, the global data and weight matrices")
}

# The sizes of the models of every sub-model, beside a row for what the global
# model holds besides them.
.global_model_size <- function(object, ...) {

  if (is.null(object[["submodels"]])) {
    stop("The global model has no sub-models yet, so there is nothing to size. ",
         "Use add_submodels() first.")
  }

  flat <- .submodel_models(object)
  models <- flat[["models"]]
  names(models) <- .submodel_model_labels(rep(names(object[["submodels"]]), flat[["sizes"]]))

  global <- as.numeric(utils::object.size(object[names(object) != "submodels"]))
  .with_global_row(bayests::expected_model_size(models, ...), global,
                   "global data and weight matrices")
}

# "US/1", "US/2", "JP/1", ...: a sub-model and the position of a model among its
# models, from the sub-model of each model in order.
.submodel_model_labels <- function(submodels) {
  position <- stats::ave(seq_along(submodels), submodels, FUN = seq_along)
  paste(submodels, position, sep = "/")
}

.with_global_row <- function(size, bytes, element) {
  global <- data.frame(model = "global", element = element, step = "",
                       draws = NA, columns = NA, bytes = bytes,
                       stringsAsFactors = FALSE)
  result <- rbind(global, as.data.frame(size))
  class(result) <- c("modelsize", "data.frame")
  result
}
