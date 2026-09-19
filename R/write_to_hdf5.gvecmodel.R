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
#' @param append logical. If \code{TRUE}, the sub-models of \code{object} are
#' written into an export that is already in \code{folder}: their files are
#' added, or replaced where the folder holds them already, and their rows are
#' merged into the manifest. The global data, the index and the weight matrices
#' of the export are left as they are, since every batch of one global model
#' describes the same ones. This is what lets a run that cannot hold all of its
#' sub-models at once estimate them in batches, writing and freeing each batch
#' before the next: the time varying models of a global model with many
#' sub-models need tens of gigabytes together and a few per batch. Defaults to
#' \code{FALSE}, which writes a fresh export.
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
#' 'gvarmodel', and so is the way to write it: a global model is estimated one
#' batch of sub-models at a time, each batch written with \code{append = TRUE}
#' and dropped before the next, which keeps in memory what one batch needs
#' rather than what the whole model does. See
#' \code{\link{write_to_hdf5.gvarmodel}} for the loop and for what reads the
#' folder afterwards.
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
                                    mc.cores = .default_cores(),
                                    append = FALSE, ...) {

  .write_gvar_to_hdf5(object = object, folder = folder, overwrite = overwrite,
                      mc.cores = mc.cores, append = append)
}
