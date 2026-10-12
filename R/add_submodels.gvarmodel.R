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
#' @param p_exogen an integer vector of the number of lags of a sub-model's
#' weakly exogenous variables, counted from the contemporaneous term onwards
#' (default is \code{p_exogen = 2}). A value of 1 uses the contemporaneous
#' variables alone, 2 adds their first lag, and 0 leaves them out altogether.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the number of lags of a sub-model's global
#' variables, counted in the same way as \code{p_exogen}. If \code{NULL}
#' (default), models do not include global variables.
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
#' @param algorithm algorithm that should be used for posterior simulation. If
#' \code{NULL} (default), the algorithm is named by \code{tvp} and
#' \code{error}. The one non-standard option is \code{"discount"}. See
#' \code{\link{create_varxsubmodel}}.
#' @param volatility_break an optional date of the data, \code{c(year, period)}
#' such as \code{c(2020, 1)}, from which every shock of the sub-model is scaled
#' by \eqn{s_t}, the volatility break of Lenza and Primiceri (2022). How large
#' the break was and how fast it faded are estimated. Available for a sub-model
#' with constant coefficients and \code{error = "wishart"} only, and refused
#' elsewhere. The date is a property of the global model rather than of one
#' country, so \code{\link{add_submodels}} is where it is normally given.
#' @param volatility_scales a whole number of at least one, how many periods
#' from \code{volatility_break} on get a scale of their own before the excess
#' decays (defaults to 3). Used only with \code{volatility_break}.
#' @param student_t logical. If \code{TRUE} the shocks of each period are
#' multivariate Student t, \eqn{u_t \sim N(0, \lambda_t \Sigma)}, so a period
#' with large shocks counts for less without anyone dating it. Available for a
#' sub-model with constant coefficients and \code{error = "wishart"} only, and
#' refused elsewhere, which includes every error correction sub-model.
#' @param delta_beta,delta_sigma numeric discount factors in \eqn{(0, 1]} of the
#' discounted model, the first governing the coefficients and the second the
#' error covariance. Both default to 1, at which the quantity they govern does
#' not move. Ignored unless \code{algorithm = "discount"}, and a vector in
#' either produces one sub-model per value. See
#' \code{\link{create_varxsubmodel}}.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' @param thin an integer thinning interval of the sampler (defaults to 1). After
#' the burn-in the sampler keeps the last of every \code{thin} draws, so it runs
#' \code{burnin + iterations * thin} draws and still keeps \code{iterations}.
#' See \code{\link{create_varxsubmodel}}.
#' @param ... further arguments passed to or from other methods.
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
#' 
#' @export
add_submodels.gvarmodel <- function(object,
                                    endogen = NULL,
                                    p_endogen = 1,
                                    exogen = NULL,
                                    p_exogen = 2,
                                    global = NULL,
                                    s = NULL,
                                    deterministic = "const",
                                    seasonal = FALSE,
                                    structural = FALSE,
                                    tvp = FALSE,
                                    error = "wishart",
                                    varsel = "none",
                                    algorithm = NULL,
                                    volatility_break = NULL,
                                    volatility_scales = 3,
                                    student_t = FALSE,
                                    delta_beta = 1,
                                    delta_sigma = 1,
                                    iterations = 10000,
                                    burnin = 2000,
                                    thin = 1,
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
                                                      algorithm = algorithm,
                                                      volatility_break = volatility_break,
                                                      volatility_scales = volatility_scales,
                                                      student_t = student_t,
                                                      delta_beta = delta_beta,
                                                      delta_sigma = delta_sigma,
                                                      iterations = iterations, burnin = burnin,
                                                      thin = thin)
  }
  
  
    
  return(object)
}
