#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list of models, usually, the output of a call to
#' \code{\link{create_gvarmodel}} in combination with \code{\link{add_submodels}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("gvar2019")
#' global_data <- gvar2019[["global_data"]]
#' submodel_data <- gvar2019[["submodel_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#' 
#' # Create global model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2013:2016)
#' 
#' object <- add_submodels(object,
#'                         p_endogen = 1,
#'                         p_exogen = 2,
#'                         global = "poil",
#'                         s = 2,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 0.0001))
#' 
#' @export
#' @method add_priors gvarmodel
add_priors.gvarmodel <- function(object, ...){
  
  # This assumes that lower level objects are of class 'modellist', which
  # is handled by bvartools-functions.
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_priors(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}