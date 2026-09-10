#' Export to HDF5 Files
#'
#' Exports the content of an object of class 'gvecmodel' to a dedicated folder.
#'
#' @param object list of class 'gvecmodel'.
#' @param folder path to the folder where the content of argument
#' \code{object} should be saved.
#' @param overwrite logical. If \code{TRUE}, an export already present in
#' \code{folder} is replaced. Defaults to \code{FALSE}, which makes the function
#' stop rather than touch it.
#' @param mc.cores the number of cores to use, i.e. at most how many sub-models
#' are written at the same time. Defaults to the number of available cores, up
#' to eight. Use \code{1} to write the sub-models one after the other.
#' Sub-models are written to separate files, so nothing is shared between the
#' workers. In contrast to the rest of the package this uses a socket cluster,
#' which also runs in parallel under Windows.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'
#' The export has the same shape as the one produced for an object of class
#' 'gvarmodel'. See \code{\link{write_to_hdf5.gvarmodel}}.
#'
#' @return The path to \code{folder}, invisibly.
#'
#' @examples
#'
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#'
#' # Set up a model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#'
#' # Export it
#' folder <- file.path(tempdir(), "bgvars-example-write-gvec")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE)
#'
#' @export
#' @method write_to_hdf5 gvecmodel
write_to_hdf5.gvecmodel <- function(object, folder, overwrite = FALSE,
                                    mc.cores = .default_cores(), ...) {

  .write_gvar_to_hdf5(object = object, folder = folder, overwrite = overwrite,
                      mc.cores = mc.cores)
}
