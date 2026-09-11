#' Printing the Data of Sub-Models
#'
#' print method for objects of class 'submodeldata'.
#'
#' @param x an object of class 'submodeldata'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The printed table, invisibly: one row per sub-model and one column
#' per variable, with a mark for each variable the sub-model provides.
#'
#' @export
print.submodeldata <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  submodel_names <- names(x)
  
  vars <- unique(unlist(lapply(x, function(y) {dimnames(y[["endogen"]])[[2]]})))
  
  availability_matrix <- matrix("", nrow = length(submodel_names), ncol = length(vars))
  dimnames(availability_matrix) <- list(submodel_names, vars)
  
  for (i in submodel_names) {
    vars_i <- dimnames(x[[i]][["endogen"]])[[2]]
    availability_matrix[i, vars_i] <- "x"
  }
  
  res <- data.frame(submodel = submodel_names, as.data.frame(availability_matrix, row.names = FALSE))
  
  print(res, digits = digits, row.names = FALSE, ...)
  invisible(res)
} 