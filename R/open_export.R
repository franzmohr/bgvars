#' Open an Exported Global Model
#'
#' Opens the export of a global model in a folder and returns a handle to it,
#' which the estimation steps of the package work through without reading the
#' sub-models into the session.
#'
#' @param folder path to a folder written by
#' \code{\link[bvartools]{write_to_hdf5}}.
#'
#' @details
#' A global model holds the posterior draws of every sub-model, and those draws
#' are what it is large in: a time varying sub-model of six variables keeps a
#' coefficient path per draw, some 0.6 GB at 3000 draws, so thirty of them need
#' tens of gigabytes at once. The export in a folder holds one file per model
#' beside a \code{model.h5} with the global data, the weight matrices and the
#' manifest, and every step of the estimation touches one model at a time.
#'
#' A handle carries the manifest and the path, not the models. The methods of
#' \code{\link{map_submodels}} -- \code{\link[bvartools]{add_priors}},
#' \code{add_initial_values}, \code{add_seed},
#' \code{add_posterior_coefficients}, \code{add_posterior_loglik},
#' \code{thin} and the two error correction transformations -- read a model,
#' change it, write it back and drop it, so what a run needs is one model per
#' worker rather than the whole model.
#'
#' The pipeline then reads as it does in memory:
#'
#' \preformatted{
#' export <- open_export(folder)
#' export <- add_priors(export, coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 1e-4))
#' export <- add_initial_values(export)
#' export <- add_seed(export, 20260919)
#' export <- add_posterior_coefficients(export, cores = 8)
#' export <- add_posterior_loglik(export, cores = 8)
#' }
#'
#' What follows the estimation reads the folder as well:
#' \code{\link{choose_best_model_from_hdf5}} ranks the models of every sub-model
#' without their coefficient draws, \code{\link{read_submodel_from_folder}}
#' reads one sub-model, and \code{\link{read_gvar_from_folder}} brings back one
#' model per sub-model to solve.
#'
#' @return A list of class 'gvarexport' with the elements \code{folder},
#' \code{manifest}, the table of \code{model.h5}, and \code{rclass}, the class
#' of the object the export was written from.
#'
#' @examples
#'
#' # Load data
#' data("gvar2023")
#' submodel_data <- select_list_elements(gvar2023[["submodel_data"]],
#'                                       c("AT", "DE", "US"))
#'
#' # Set up a model and write it
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
#' folder <- file.path(tempdir(), "bgvars-example-export")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
#'
#' # Estimate it through the folder
#' export <- open_export(folder)
#' export <- add_priors(export, coef = list(v_i = 1),
#'                      sigma = list(df = 3, scale = 1e-8))
#' export <- add_initial_values(export)
#' export <- add_posterior_coefficients(export)
#' export
#'
#' @family export
#' @export
open_export <- function(folder) {

  if (!dir.exists(folder)) {
    stop("Folder ", folder, " does not exist.")
  }
  filename <- file.path(folder, "model.h5")
  if (!file.exists(filename)) {
    stop("Folder ", folder, " holds no model.h5, so it is not an export of a ",
         "global model.")
  }

  file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE)

  rclass <- if (file$attr_exists("rclass")) hdf5r::h5attr(file, "rclass") else NA_character_
  manifest <- if (file$exists("submodels")) {
    as.data.frame(file[["submodels"]]$read(), stringsAsFactors = FALSE)
  } else {
    NULL
  }
  file$close_all()

  if (is.null(manifest) || nrow(manifest) == 0) {
    stop("The export in ", folder, " holds no sub-models.")
  }
  missing <- !file.exists(file.path(folder, manifest[["file"]]))
  if (any(missing)) {
    stop(sum(missing), " of the ", nrow(manifest), " files the manifest names ",
         "are not in ", folder, ", so the export is incomplete.")
  }

  structure(list(folder = folder, manifest = manifest, rclass = rclass),
            class = c("gvarexport", "list"))
}

#' @rdname open_export
#' @param x an object of class 'gvarexport'.
#' @param ... further arguments passed to or from other methods.
#' @export
print.gvarexport <- function(x, ...) {

  manifest <- x[["manifest"]]
  submodels <- unique(manifest[["submodel"]])
  cat("Export of a global model in", x[["folder"]], "\n")
  cat(nrow(manifest), "models of", length(submodels), "sub-models\n")

  specification <- unique(data.frame(
    type = manifest[["type"]], tvp = manifest[["tvp"]],
    error = manifest[["error"]], varsel = manifest[["varsel"]],
    iterations = manifest[["iterations"]], burnin = manifest[["burnin"]],
    stringsAsFactors = FALSE))
  cat("\n")
  print(specification, row.names = FALSE)

  invisible(x)
}

# The manifest of an export, written back to model.h5 after the models it
# describes have changed on disk.
.write_export_manifest <- function(object) {

  file <- hdf5r::h5file(file.path(object[["folder"]], "model.h5"), mode = "a")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE)
  if (file$exists("submodels")) {
    file$link_delete("submodels")
  }
  file[["submodels"]] <- object[["manifest"]]
  file$close_all()

  invisible(object)
}
