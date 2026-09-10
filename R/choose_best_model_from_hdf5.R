#' Choose Best Sub-Model
#' 
#' Obtains model selection criteria for each model and gives
#' a table of the best performing model for each sub-model.
#' 
#' @param folder path to a folder holding an export, that is, a \code{model.h5}
#' and a \code{submodels} directory beside it.
#' @param criterion selection criterion used to choose the best model using
#' \code{\link[bvartools]{selection_criteria}}.
#' 
#' @details
#' Only the parts of a model file that enter a selection criterion are read,
#' that is, the draws of the log-likelihood, the forecast errors and the
#' handful of model specifications the criteria are calculated from. The
#' posterior draws of the coefficients make up most of an exported model and
#' are never touched, so the function does not need to hold a whole sub-model
#' in memory to rank it.
#' 
#' @return A data frame containing the postion of the best performing model
#' for each sub-model.
#' 
#' @export
choose_best_model_from_hdf5 <- function(folder, criterion) {
  
  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " does not exist."))
  }
  
  filename_model <- file.path(folder, "model.h5")
  if (!file.exists(filename_model)) {
    stop(paste0("Folder ", folder, " does not contain a model.h5."))
  }
  
  model_file <- hdf5r::h5file(filename_model, mode = "r")
  on.exit(if (model_file$is_valid) model_file$close_all(), add = TRUE)
  
  # Sub-models ----
  #
  # Through the manifest rather than through the directory tree: which
  # sub-model a file belongs to is recorded there, and a folder that has been
  # moved or renamed since must not change the answer.
  manifest <- NULL
  if ("submodels" %in% names(model_file)) {
    manifest <- model_file[["submodels"]]$read()
  }
  model_file$close_all()
  
  if (is.null(manifest)) {
    stop("No submodels available.")
  }
  
  submodel <- unique(manifest[, "submodel"])
  
  result <- data.frame("submodel" = submodel, "position" = NA_integer_)
  for (i in 1:length(submodel)) {
    
    files <- manifest[manifest[, "submodel"] == submodel[i], "file"]
    
    # The criteria are obtained file by file, so that at no point more than one
    # model's draws of the log-likelihood are in memory.
    sc <- vector(mode = "list", length = length(files))
    for (j in seq_along(files)) {
      sc[[j]] <- bvartools::selection_criteria(
        .read_selection_data_from_hdf5(file.path(folder, files[j])))
    }
    class(sc) <- c("selcritlist", "list")
    
    result[i, "position"] <- bvartools::choose_best_model(sc, criterion = criterion)
  }
  
  return(result)
}

#' Read the Inputs of a Selection Criterion from an HDF5 File
#' 
#' Reads the parts of an exported model that
#' \code{\link[bvartools]{selection_criteria}} works on and returns them as a
#' model object of the class the file was written with.
#' 
#' @param filename path to a single exported model.
#' 
#' @details
#' The training data only enters the criteria through its dimensions, so its
#' shape is taken from the file and the data itself is not read. The same holds
#' for the specification of the model, of which only the entries a criterion is
#' calculated from are read.
#' 
#' Handles are closed one by one instead of through \code{close_all}, which
#' walks every handle the file has ever opened and, over the hundreds of files
#' a global model consists of, costs more time than reading the draws.
#' 
#' @return A list of the class recorded in the file.
#' 
#' @noRd
.read_selection_data_from_hdf5 <- function(filename) {
  
  if (!file.exists(filename)) {
    stop(paste0("Model file ", filename, " does not exist."))
  }
  
  h5 <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5$is_valid) h5$close_all(), add = TRUE)
  handles <- list()
  
  # Specification ----
  #
  # Everything else the group records describes the model rather than entering
  # a criterion, and is available from the manifest for the models it is
  # wanted for.
  group <- h5[["model"]]
  handles <- c(handles, group)
  used_attrs <- c("k", "p", "m", "s", "n", "rank", "h", "structural", "rclass")
  model <- list()
  for (attr_name in intersect(used_attrs, hdf5r::h5attr_names(group))) {
    model[[attr_name]] <- hdf5r::h5attr(group, attr_name)
  }
  
  object <- list("model" = model)
  
  # Training data ----
  #
  # Only the number of observations and the number of regressors are used.
  if (h5$exists("data/train/y")) {
    dataset <- h5[["data/train/y"]]
    handles <- c(handles, dataset)
    object[["data"]][["train"]][["y"]] <- matrix(NA_real_, dataset$dims[1], 1L)
  }
  if (h5$exists("data/train/x")) {
    dataset <- h5[["data/train/x"]]
    handles <- c(handles, dataset)
    object[["data"]][["train"]][["x"]] <- matrix(NA_real_, 1L, dataset$dims[2])
  }
  
  # Posterior ----
  if (h5$exists("posterior/loglik")) {
    dataset <- h5[["posterior/loglik"]]
    handles <- c(handles, dataset)
    object[["posterior"]][["loglik"]] <- dataset$read()
  }
  if (h5$exists("posterior/forecast_errors")) {
    dataset <- h5[["posterior/forecast_errors"]]
    handles <- c(handles, dataset)
    object[["posterior"]][["forecast_errors"]] <- dataset$read()
    # Forecast errors are the only output that is labelled by variable.
    endogen <- hdf5r::h5attr(group, "endogen")
    object[["data"]][["original"]][["endogen"]] <-
      matrix(NA_real_, 1L, length(endogen), dimnames = list(NULL, endogen))
  }
  
  for (handle in handles) {
    handle$close()
  }
  h5$close()
  
  class(object) <- model[["rclass"]]
  
  return(object)
}
