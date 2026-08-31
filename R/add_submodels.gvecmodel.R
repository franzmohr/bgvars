#' Add Sub-Models to a Global Model
#'
#' Generate and add VECX sub-models for a global vector error correction model.
#'
#' @param object a list of class 'gvecmodel'.
#' @param endogen character vector of variables that should enter each sub-model
#' as endogenous variables, if they are available for the respective sub-model.
#' @param p_endogen an integer vector of the lag order of endogenous variables
#' in the (levels) VAR. Thus, the resulting model's lag will be \code{p_endogen - 1}.
#' @param exogen character vector of variables that should enter each sub-model
#' as exogenous variables.
#' @param p_exogen an integer vector of the lag order of a sub-model's weakly
#' exogenous variables in the (levels) VAR. Thus, the resulting model's lag will
#' be \eqn{p_exogen - 1}.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the lag order of a sub-model's global variables
#' in the (levels) VAR. Thus, the resulting model's lag will be \code{s - 1}.
#' If \code{NULL} (default), models do not include global variables.
#' @param r an integer vector of the cointegration rank.
#' @param const a character specifying whether a constant term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param trend a character specifying whether a trend term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param seasonal a character specifying whether seasonal dummies should be included in the error
#' correction term (\code{"restricted"}) or in the non-cointegreation term as \code{"unrestricted"}
#' variables. If \code{NULL} (default) no seasonal terms will be added. The amount of dummy variables
#' will be automatically detected and depends on the frequency of the time-series object provided
#' in \code{data}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model. Defaults to \code{FALSE}.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param algorithm algorithm that should be used for posterior simulation. If \code{NULL}
#' (default), standard algorithms will be used. See 'Details' for available
#' non-standard options.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details
#' Available specifications for argument \code{algorithm} are:
#' \itemize{
#'  \item{\code{"KLGS2010"}: Algorithm proposed in Koop, León-González & Strachan (2010).}
#' }
#' 
#' 
#' @examples
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
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
#' object <- add_submodels(object, r = 1)
#' 
#' 
#' 
#' @export
add_submodels.gvecmodel <- function(object,
                                    endogen = NULL,
                                    p_endogen = 1,
                                    exogen = NULL,
                                    p_exogen = 1,
                                    global = NULL,
                                    s = NULL,
                                    r = NULL,
                                    const = NULL,
                                    trend = NULL,
                                    seasonal = NULL,
                                    structural = FALSE,
                                    error = "wishart",
                                    tvp = FALSE,
                                    varsel = "none",
                                    algorithm = NULL,
                                    iterations = 10000,
                                    burnin = 2000,
                                    ...){
  
  # submodel <- "BR"; endogen = c("y","p", "rs"); p_endogen = 1; exogen = c("y", "p"); p_exogen = 1; global = "poil";
  # s = 0; deterministic = "const"; seasonal = FALSE;  structural = FALSE;
  # tvp = FALSE; error = "wishart"; varsel = "none"; iterations = 10000; burnin = 2000
  # for (i in submodels) {
  #   temp <- create_varxsubmodel(i,
  #   object = object,
  #   endogen = endogen, p_endogen = p_endogen,
  #   exogen = exogen, p_exogen = p_exogen,
  #   global = global, s = s,
  #   deterministic = deterministic,
  #   seasonal = seasonal,
  #   structural = structural,
  #   tvp = tvp, error = error, varsel = varsel,
  #   iterations = iterations, burnin = burnin)
  # }
  
  index <- object[["global"]][["index"]]
  submodels <- unique(index[, "submodel"])
  
  for (i in submodels) {
    object[["submodels"]][[i]] <- create_vecxsubmodel(object, submodel = i,
                                                      endogen = endogen, p_endogen = p_endogen,
                                                      exogen = exogen, p_exogen = p_exogen,
                                                      global = global, s = s,
                                                      r = r,
                                                      const = const, trend = trend, seasonal = seasonal,
                                                      structural = structural,
                                                      error = error, tvp = tvp, varsel = varsel,
                                                      algorithm = algorithm,
                                                      iterations = iterations, burnin = burnin)
  }
  
  
    
  return(object)
}
