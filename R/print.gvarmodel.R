#' Printing Model Information
#'
#' print method for objects of class 'gvarmodel'.
#'
#' @param x an object of class 'gvarmodel'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The object that was printed, invisibly. A table of the variables of
#' every sub-model is printed as a side effect, with a mark for each variable
#' that the sub-model contains.
#'
#' @export
print.gvarmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  # Used variables
  avail_vars <- x[["global"]][["index"]][, c("submodel", "variable")]
  
  submodels <- unique(avail_vars[, "submodel"])
  vars <- unique(avail_vars[, "variable"])
  
  availability_matrix <- matrix("", nrow = length(submodels), ncol = length(vars))
  dimnames(availability_matrix) <- list(submodels, vars)
  
  for (i in submodels) {
    vars_i <- avail_vars[avail_vars[, "submodel"] == i, "variable"]
    availability_matrix[i, vars_i] <- "x"
  }
  
  res <- data.frame(submodel = submodels, as.data.frame(availability_matrix, row.names = FALSE))
  
  cat("\nGlovar Vector Autoregressive Model\n\n")
  
  cat("\nVariables in sub-models:\n\n")
  print(res, digits = digits, row.names = FALSE, ...)
  
  if (!is.null(x[["global"]][["exogen"]])) {
    res <- dimnames(x[["global"]][["exogen"]])[[2]]
    cat("\nGlobal variables:\n\n")
    cat("\t", paste0(res, collapse = "\n\t"), sep = "")
  }
  
  # Number of models per sub-model
  n_models <- lapply(x[["submodels"]], length)
  
  if (length(n_models) > 0) {
    stop("Update print.gvarmodel to handle cases, where submodels are available in the model.")
  }
  
} 