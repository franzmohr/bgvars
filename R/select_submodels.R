#' Submodel Selection
#' 
#' Creates a list of submodels of a GVAR model, which are used to solve the
#' global model.
#' 
#' @param object a list containing the posterior draws of the submodels of a
#' GVAR model. Usually, the result of a call to \code{\link[bvartools]{draw_posterior}}.
#' @param models a vector of integers with the position of the submodels, which
#' should be used in the global model.
#' 
#' @export
select_submodels <- function(object, models) {
  
  obj_class <- class(object)
  if (!"submodelestlist" %in% obj_class) {
    stop("Argument 'object' must be of class 'submodelestlist'.")
  }
  
  obj_names <- unique(names(object))
  if (!all(names(object)[models] %in% obj_names)) {
    stop("Selection in argument 'models' would lead to an incomplete global model.")
  }
  if (any(table(names(object)[models]) > 1)) {
    stop("Selection in argument 'models' would lead to more than one entity or country being used in the global model.")
  }

  result <- list()
  for (i in 1:length(models)) {
    result[[i]] <- object[[models[i]]]
  }
  names(result) <- names(object)[models]
  
  class(result) <- obj_class
  
  return(result)
}