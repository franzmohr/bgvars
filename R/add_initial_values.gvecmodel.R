#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a list of models by passing each element to
#' the respective method.
#'
#' @param object an object of class 'gvecmodel', usually, the output of a call to
#' \code{\link{create_gvecmodel}} in combination with \code{\link{add_submodels}}
#' and \code{\link[bvartools]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvecmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("EA", "US"))
#' 
#' # Create empty model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#' 
#' # Create sub-models
#' object <- add_submodels(object,
#'                         p_endogen = 1, p_exogen = 1,
#'                         global = "poil",
#'                         s = 1,
#'                         r = 1,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#'                      
#' # Add initial values
#' object <- add_initial_values(object)
#' 
#' @export
#' @method add_initial_values gvecmodel
add_initial_values.gvecmodel <- function(object, ...){
  
  # This assumes that lower level objects are of class 'modellist', which
  # is handled by bvartools-functions.
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_initial_values(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}