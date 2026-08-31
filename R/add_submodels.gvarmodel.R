#' Add Sub-Models to a Global Model
#'
#' Generate and add VARX sub-models for a global vector autoregressive model.
#'
#' @param object a list of class 'gvarmodel'.
#' @param endogen character vector of variables that should enter each sub-model
#' as endogenous variables, if they are available for the respective sub-model.
#' @param p_endogen an integer vector of the lag order (default is \code{p_endogen = 1})
#' of a sub-model's endogenous variables.
#' @param exogen character vector of variables that should enter each sub-model
#' as weakly exogenous variables.
#' @param p_exogen an integer vector of the lag order (default is \code{p_exogen = 1})
#' of a sub-model's weakly exogenous variables.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the lag order of a sub-model's global variables.
#' If \code{NULL} (default), models do not include global variables.
#' @param deterministic a character specifying which deterministic terms should
#' be included. Available values are \code{"none"}, \code{"const"} (default) for an intercept,
#' \code{"trend"} for a linear trend, and \code{"both"} for an intercept with a linear trend.
#' @param seasonal logical. If \code{TRUE}, seasonal dummy variables are
#' generated as additional deterministic terms. The amount of dummies depends on the frequency of the
#' time-series object provided in \code{object}. Defaults to \code{FALSE}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model. Defaults to \code{FALSE}.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @examples
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]]
#' global_data <- gvar2023[["global_data"]]
#' 
#' # Create empty model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#' 
#' # Create sub-models
#' object <- add_submodels(object,
#'                         endogen = c("y","p", "rs"), p_endogen = 1,
#'                         exogen = c("y", "p"), p_exogen = 1,
#'                         global = "poil", s = 0,
#'                         deterministic = "const", seasonal = FALSE,
#'                         structural = FALSE, tvp = FALSE,
#'                         error = "wishart", varsel = "none",
#'                         iterations = 10000, burnin = 2000)
#' 
#' 
#' 
#' @export
add_submodels.gvarmodel <- function(object,
                                    endogen = NULL,
                                    p_endogen = 1,
                                    exogen = NULL,
                                    p_exogen = 1,
                                    global = NULL,
                                    s = NULL,
                                    deterministic = "const",
                                    seasonal = FALSE,
                                    structural = FALSE,
                                    tvp = FALSE,
                                    error = "wishart",
                                    varsel = "none",
                                    iterations = 10000,
                                    burnin = 2000,
                                    ...){
  
  submodels <- unique(object[["global"]][["index"]][, "submodel"])
  
  for (i in submodels) {
    object[["submodels"]][[i]] <- create_varxsubmodel(object, submodel = i,
                                                      endogen = endogen, p_endogen = p_endogen,
                                                      exogen = exogen, p_exogen = p_exogen,
                                                      global = global, s = s,
                                                      deterministic = deterministic,
                                                      seasonal = seasonal,
                                                      structural = structural,
                                                      tvp = tvp, error = error, varsel = varsel,
                                                      iterations = iterations, burnin = burnin)
  }
  
  
    
  return(object)
}
