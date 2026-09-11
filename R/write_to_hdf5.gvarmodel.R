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
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#'
#' # Set up a model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#'                            
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2014:2016)
#'                               
#' # Create sub-models
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp", "r"),
#'                         p_endogen = 1,
#'                         exogen = c("y", "Dp", "r"),
#'                         p_exogen = 2,
#'                         global = "poil",
#'                         s = 2,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 0.0001))
#'                      
#' # Add initial values
#' object <- add_initial_values(object)
#'
#' # Export models
#' folder <- file.path(tempdir(), "bgvars-example-write-gvar")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE)
#'
#' @export
#' @method write_to_hdf5 gvarmodel
write_to_hdf5.gvarmodel <- function(object, folder, overwrite = FALSE,
                                    mc.cores = .default_cores(), ...) {

  .write_gvar_to_hdf5(object = object, folder = folder, overwrite = overwrite,
                      mc.cores = mc.cores)
}
