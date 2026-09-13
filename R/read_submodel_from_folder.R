#' Import Sub-Models from a Folder of HDF5 Files
#' 
#' Imports every model stored below a sub-model's folder.
#' 
#' @param folder Path to the root folder of a global model. Usually, the same
#' path used as argument \code{"folder"} during a call to
#' \code{\link{write_to_hdf5.gvarmodel}}.
#' @param submodel character of the sub-model that should be imported.
#' 
#' @details
#' The function allows to import the data from a single sub-model of a global
#' model. It's main purpose is to avoid reading all models of a global model at
#' once in order to save memory.
#' 
#' @return A list of class 'modellist'.
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
#' folder <- file.path(tempdir(), "bgvars-example-read-submodel")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE)
#' 
#' # Import
#' submodel <- read_submodel_from_folder(folder = folder, submodel = "AT")
#' 
#' 
#' @export
read_submodel_from_folder <- function(folder, submodel) {
  
  if (!dir.exists(folder)) {
    stop("Specified folder does not exist.")
  }
  
  if (length(submodel) > 1) {
    stop("Argument 'submodel' must have length 1.")
  }
  
  if (!"character" %in% class(submodel)) {
    stop("Argument 'submodel' must be of class 'character'.")
  }
  
  # Through the manifest where the export has one, as read_gvar_from_folder
  # does. The position of a model -- which is what choose_best_model_from_hdf5
  # returns -- is its row among the rows of its sub-model in the manifest.
  # Walking the directory instead orders the files by name as text, which
  # agrees with the manifest only while every name has the same number of
  # digits.
  filename_model <- file.path(folder, "model.h5")
  if (file.exists(filename_model)) {
    model_file <- hdf5r::h5file(filename_model, mode = "r")
    on.exit(if (model_file$is_valid) model_file$close_all(), add = TRUE)
    manifest <- NULL
    if ("submodels" %in% names(model_file)) {
      manifest <- model_file[["submodels"]]$read()
    }
    model_file$close_all()

    if (!is.null(manifest)) {
      rows <- manifest[manifest[, "submodel"] == submodel, , drop = FALSE]
      if (nrow(rows) == 0) {
        stop("Submodel ", submodel, " is not available. The export does not contain it.")
      }

      result <- .read_submodels(folder, rows)[[submodel]]
      # Named after the files, as bvartools::read_models_from_folder names them.
      names(result) <- sub("[.]h5$", "", basename(rows[, "file"]))
      return(result)
    }
  }

  path_to_submodel <- file.path(folder, "submodels", submodel)

  if (!dir.exists(path_to_submodel)) {
    stop("Submodel ", submodel, " is not available. Its folder does not exist.")
  }

  return(bvartools::read_models_from_folder(path_to_submodel))

}