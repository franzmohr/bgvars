#' Export to HDF5 Files
#'
#' Exports the content of an object of class 'gvarmodel' to a dedicated folder.
#'
#' @param object list of class 'gvarmodel'.
#' @param folder path to the folder where the content of argument
#' \code{object} should be saved.
#' @param overwrite logical. If \code{TRUE}, an export already present in
#' \code{folder} is replaced. Defaults to \code{FALSE}, which makes the function
#' stop rather than touch it.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'
#' The export consists of one file for the global model and one file per
#' sub-model:
#'
#' \describe{
#'   \item{\code{model.h5}}{The global data in group \code{global}, one weight
#'   matrix per sub-model in group \code{weights}, and a \code{submodels} table
#'   describing every sub-model file.}
#'   \item{\code{submodels/<sub-model>/<nnn>.h5}}{One estimable model each,
#'   numbered in the order in which \code{\link{add_submodels}} produced them.
#'   The \code{submodels} table in \code{model.h5} says which specification each
#'   of them holds.}
#' }
#'
#' Sub-models are kept in separate files because that is the unit an external
#' sampler can work on independently: HDF5 does not support concurrent writers
#' to one file.
#'
#' @return The path to \code{folder}, invisibly.
#'
#' @examples
#'
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]]
#' global_data <- gvar2023[["global_data"]]
#'
#' # Set up a model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#'
#' # Export it
#' folder <- file.path(tempdir(), "gvar")
#' dir.create(folder)
#' write_to_hdf5(object, folder = folder)
#'
#' @export
#' @method write_to_hdf5 gvarmodel
write_to_hdf5.gvarmodel <- function(object, folder, overwrite = FALSE, ...) {

  .write_gvar_to_hdf5(object = object, folder = folder, overwrite = overwrite)
}
