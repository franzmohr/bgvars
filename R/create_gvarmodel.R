#' Create a Global Vector Autoregressive Model
#'
#' Initialises an object of class 'gvarmodel', which contains all elements necessary
#' to set up, estimate and evaluate a Bayesian GVAR model.
#'
#' @param submodel_data an object of class 'submodeldata', which is a named list
#' of time-series objects of with endogenous data for sub-models.
#' @param global_data a named time-series object of global data.
#' 
#' @details
#' The function creates a list of class 'gvarmodel', which provides the basic
#' structure for the creation of sub-models, prior specification, initial value
#' generation, posterior simulation, model evaluation and structural analysis.
#' 
#' @details
#' The deterministic terms of the global model --- a constant, a linear trend
#' and, for data of a frequency above one, a set of seasonal dummies --- are
#' built here, on the time axis of the data, and stored as element
#' \code{deterministic}. Every sub-model takes the ones it uses from this
#' series rather than building its own.
#'
#' The reason is the trend. A sub-model that counted a trend from its own first
#' observation would not use the same regressor as its neighbours: a sub-model
#' with fewer lags keeps more of the early observations, so its trend would be
#' shifted against that of a sub-model with more lags. Each of them would be
#' internally consistent, since a shift of the trend is absorbed by the
#' constant, but the global model has a single trend regressor, and
#' \code{\link{submodels_to_gvar}} could not stack sub-models that disagree on
#' what it is.
#'
#' @return A list of class 'gvarmodel', which consists of the following elements:
#' \describe{
#'   \item{\strong{global}}{A named list containing the data objects of the global
#'   model: the endogenous variables of all units, the global variables, the
#'   deterministic terms and an index of the variables. See 'Details'.}
#'   \item{\strong{weights}}{The list entry, where the submodel weight matrices will go.
#'   This entry is empty after the execution of \code{create_gvarmodel}, but will be
#'   updated with \code{\link{add_weight_matrices}}.} 
#'   \item{\strong{submodels}}{The list entry, where the estimated submodels will go.
#'   This entry is empty after the execution of \code{create_gvarmodel}, but will be
#'   updated with \code{\link{add_submodels}}.}
#' }
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
#' @export
create_gvarmodel <- function(submodel_data, global_data = NULL){
  
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
                                 "deterministic" = .global_deterministic(endogen),
                                 "index" = index),
                 "weights" = NULL,
                 "submodels" = NULL)
  
  class(result) <- list("gvarmodel", "list")  
  
  return(result)
}
