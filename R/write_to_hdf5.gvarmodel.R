#' Export to HDF5 File
#'
#' Exports the content of an object of class 'gvarmodel' to a dedicated folder.
#'
#' @param object list of class 'gvarmodel'.
#' @param folder path to the folder where the content of argument
#' \code{object} should be saved.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' 
#' @export
#' @method write_to_hdf5 gvarmodel
write_to_hdf5.gvarmodel <- function(object, folder, ...) {
  
  if (!dir.exists(folder)) {
    stop(paste("Folder", folder, "does not exist."))
  }
  
  try({
    
    # **************************************************************************
    # Save global data ----
    filename_global <- file.path(folder, "global.h5")
    
    if (file.exists(filename_global)) {
      stop(paste0("File ", filename_global, " already exists."))
    }
    
    global <- hdf5r::h5file(filename_global, mode = "a")
    for (i in c("endogen", "exogen")) {
      if (!is.null(object[["global"]][[i]])) {
        global[[i]] <- object[["global"]][[i]]
        hdf5r::h5attr(global[[i]], "variables") <- dimnames(object[["global"]][[i]])[[2]]
        hdf5r::h5attr(global[[i]], "tsp") <- stats::tsp(object[["global"]][[i]])
      } 
    }
    global[["index"]] <- object[["global"]][["index"]]
    
    # Close file
    global$close_all()
    
    # **************************************************************************
    # Save weight data ----
    filename_weights <- file.path(folder, "weights.h5")
    
    if (file.exists(filename_weights)) {
      stop(paste0("File ", filename_weights, " already exists."))
    }
    
    weights <- hdf5r::h5file(filename_weights, mode = "a")
    for (i in names(object[["weights"]])) {
      weights[[i]] <- object[["weights"]][[i]]
    }
    
    # Close file
    weights$close_all()
    
    # **************************************************************************
    # Save sub-model data ----
    
    # Create submodel folder
    path_submodels <- file.path(folder, "submodels")
    if (!dir.exists(path_submodels)) {
      dir.create(path_submodels) 
    }
    for (s_i in names(object[["submodels"]])) {
      path_submodel_i <- file.path(path_submodels, s_i)
      if (!dir.exists(path_submodel_i)) {
        dir.create(path_submodel_i) 
      }
      bvartools::write_to_hdf5(object = object[["submodels"]][[s_i]], folder = path_submodel_i)
    }
  })
}