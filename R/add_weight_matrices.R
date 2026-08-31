#' Add Weight Matrices
#' 
#' Adds weight matrices for each sub-model in an object of class 'gvarmodel'
#' or 'gvecmodel'.
#' 
#' @param object an object of class 'gvarmodel' or 'gvecmodel'. Usually, the output of a call
#' to \code{\link{create_gvarmodel}} or \code{\link{create_gvecmodel}}.
#' @param submodel_data a list of class 'submodeldata' consisting of two
#' time-series objects, \code{endogen} and \code{weights}, which contain
#' sub-model-specific observations of endogenous variables and weights, respectively.
#' @param period either an integer for time varying weights or a numeric
#' vector specifying the periods in \code{weights} that
#' should be used to calculate constant weights. See 'Details'.
#' 
#' @details 
#' 
#' The function creates of sub-model-specific weight matrices.
#' If a numeric vector is provided as \code{period}, the function calculates
#' weights based on the sums over the specified periods. If an integer
#' is proved, the weights are constructed from rolling sums over the last
#' \code{period} periods. If a sub-model series begins earlier than its
#' corresponding weight series, the sums over the first \code{period}
#' observations of a sub-model's weight data are used until the periods match.
#' 
#' @return A list class 'gvarmodel' or 'gvecmodel'.
#' 
#' @examples
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]] # Country series including weights
#' global_data <- gvar2023[["global_data"]] # Global commodities data
#' 
#' # Create 'gvarmodel' object
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Generate weight matrices as 3 year, rolling window averages
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#' 
#' 
#' @export
add_weight_matrices <- function(object, submodel_data, period){
  
  # Check input
  .check_submodeldata(submodel_data)
  
  index <- object[["global"]][["index"]]
  submodel_names <- unique(index[, "submodel"])
  
  weights <- lapply(submodel_data, .create_weights_submodel, period = period)
  
  # Ensure that the order of weights is consistent with order of sub-models themselves
  # This ensures that weights are attributed correctly later.
  for (s_i in submodel_names) {
    weights[[s_i]] <- weights[[s_i]][, submodel_names]
  }
  
  # Add final weight matrices
  object[["weights"]] <- lapply(submodel_names, .create_weight_matrices, index = index, weights = weights)
  names(object[["weights"]]) <- submodel_names
  
  return(object)
}
