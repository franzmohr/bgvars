#' Create a Global Vector Error Correction Model
#'
#' Initialises an object of class 'gvecmodel', which contains all elements necessary
#' to set up, estimate and evaluate a Bayesian GVEC model.
#'
#' @param submodel_data an object of class 'submodeldata', which is a named list
#' of time-series objects of with endogenous data for sub-models.
#' @param global_data a named time-series object of global data.
#' 
#' @details
#' The function creates a list of class 'gvecmodel', which provides the basic
#' structure for the creation of sub-models, prior specification, initial value
#' generation, posterior simulation, model evaluation and structural analysis.
#' 
#' @return A list of class 'gvecmodel', which consists of the following elements:
#' \describe{
#'   \item{\strong{global}}{A named list containing the data objects of the global model.}
#'   \item{\strong{weights}}{The list entry, where the sub-model weight matrices will go.
#'   This entry is empty after the execution of \code{create_gvecmodel}, but will be
#'   updated with \code{\link{add_weight_matrices}}.} 
#'   \item{\strong{submodels}}{The list entry, where the estimated sub-models will go.
#'   This entry is empty after the execution of \code{create_gvecmodel}, but will be
#'   updated with \code{\link{add_submodels}}.}
#' }
#' 
#' @examples
#' 
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]] # Country series including weights
#' global_data <- dees2007[["global_data"]] # Global commodities data
#' 
#' # Create 'gvarmodel' object
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#'                            
#' @export
create_gvecmodel <- function(submodel_data, global_data = NULL){
  
  if (!"submodeldata" %in% class(submodel_data)) {
    stop("Argument 'submodel_data' must be of class 'submodeldata'.")
  }
  
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
      dimnames(global_data)[[2]] <- "global"
    }
  }
  
  result <- list("global" = list("endogen" = endogen,
                                 "exogen" = global_data,
                                 "index" = index),
                 "weights" = NULL,
                 "submodels" = NULL)
  
  class(result) <- list("gvecmodel", "list")  
  
  return(result)
}
