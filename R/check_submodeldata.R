

.check_submodeldata <- function(submodel_data) {
  
  submodel_names <- names(submodel_data)
  if (is.null(submodel_data)) {
    stop("Elements in object submodel_data are not named.")
  }
  
  for (s_i in submodel_names) {
    if (is.null(submodel_data[[s_i]][["endogen"]])) {
      stop(paste0("Sub-model ", s_i, " does not contain element 'endogen'."))
    }
    if (!"ts" %in% class(submodel_data[[s_i]][["endogen"]])) {
      stop(paste0("Element 'endogen' in sub-model ", s_i, " must be a time-series object."))
    }
    if (is.null(submodel_data[[s_i]][["weights"]])) {
      stop(paste0("Sub-model ", s_i, " does not contain element 'weights'."))
    }
    if (!"ts" %in% class(submodel_data[[s_i]][["weights"]])) {
      stop(paste0("Element 'weights' in sub-model ", s_i, " must be a time-series object."))
    }
    if (!all(submodel_data[[s_i]][["weights"]][, s_i] == 0)) {
      stop(paste0("Within element 'weights' of sub-model ", s_i, " a sub-model's own weights must be zero."))
    }
  }
}