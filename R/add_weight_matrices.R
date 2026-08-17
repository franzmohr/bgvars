#' Add Weight Matrices
#'  
#' Generic function used to generate weight matrices and add them to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_weight_matrices <- function (object, ...) {
 UseMethod("add_weight_matrices")
}