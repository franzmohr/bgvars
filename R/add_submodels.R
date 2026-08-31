#' Add Sub-Models to Global Model
#'  
#' Generic function used to generate sub-models and add them to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_submodels <- function (object, ...) {
 UseMethod("add_submodels")
}