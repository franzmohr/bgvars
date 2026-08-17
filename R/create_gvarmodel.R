#' Create a Global Vector Autoregressive Model
#'
#' Initialises an object of class 'gvarmodel', which contains all elements necessary
#' to set up, estimate and evaluate a Bayesian GVAR model.
#'
#' @param submodel_data a named list of time-series objects of country-specific data.
#' @param global_data a named time-series object of global data.
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]] # Country series including weights
#' global_data <- gvar2023[["global_data"]] # Global commodities data
#' 
#' # Create 'gvarmodel' object
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' 
#' 
#' 
#' @export
create_gvarmodel <- function(submodel_data, global_data = NULL){
  
  submodels <- names(submodel_data)
  
  # Collect submodel data for global model
  endogen <- NULL
  index <- NULL
  for (s_i in submodels) {
    endogen_i <- submodel_data[[s_i]][["endogen"]]
    endogen <- cbind(endogen, endogen_i)
    
    index_i <- data.frame(submodel = s_i, variable = dimnames(endogen_i)[[2]])
    index <- rbind(index, index_i)
    
    rm(list = c("endogen_i", "index_i"))
  }
  index[, "index"] <- apply(index, 1, function(x) {paste0(x["submodel"], "_", x["variable"])})
  index[, "id"] <- 1:nrow(index)
  dimnames(endogen)[[2]] <- index[, "index"]
  
  # Check global model
  if (!is.null(global_data)) {
    if (!"ts" %in% class(global_data)) {
      stop("Argument 'global_data' must be of class 'ts'.")
    }
    
    if (is.null(dimnames(global_data))) {
      tsp_global <- stats::tsp(global_data)
      # If 'global_data' is a simple ts object, transform it into a matrix object
      # to keep variable name information
      global_data <- stats::ts(as.matrix(global_data), class = c("mts", "ts", "matrix"))
      stats::tsp(global_data) <- tsp_global
    }
    dimnames(global_data)[[2]] <- "global"
  }
  
  result <- list("global" = list("endogen" = endogen,
                                 "exogen" = global_data,
                                 "index" = index),
                 "weights" = NULL,
                 "submodels" = NULL)
  
  class(result) <- list("gvarmodel", "list")  
  
  return(result)
}
