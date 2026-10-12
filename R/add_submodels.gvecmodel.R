#' Add Sub-Models to a Global Model
#'
#' Generate and add VECX sub-models for a global vector error correction model.
#'
#' @param object a list of class 'gvecmodel'.
#' @param endogen character vector of variables that should enter each sub-model
#' as endogenous variables, if they are available for the respective sub-model.
#' @param p_endogen an integer vector of the lag order of the endogenous
#' variables in the corresponding VAR model in levels (default is
#' \code{p_endogen = 1}). The sub-models hold \code{p_endogen - 1} lagged
#' differences of them. See \code{\link{create_vecxsubmodel}}.
#' @param exogen character vector of variables that should enter each sub-model
#' as exogenous variables.
#' @param p_exogen an integer vector of the lag order of the weakly exogenous
#' variables in the corresponding VARX sub-model in levels, counted from the
#' contemporaneous term onwards as in \code{\link{add_submodels.gvarmodel}}
#' (default is \code{p_exogen = 2}). The sub-models hold \code{p_exogen - 1}
#' differences of them beside their lagged levels in the error correction term.
#' See \code{\link{create_vecxsubmodel}}.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the lag order of the global variables in the
#' corresponding VARX sub-model in levels, counted in the same way as
#' \code{p_exogen}. If \code{NULL} (default), models do not include global
#' variables.
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
#' @param delta_beta,delta_sigma numeric discount factors in \eqn{(0, 1]} of the
#' discounted model, the first governing the coefficients and the second the
#' error covariance. Both default to 1, at which the quantity they govern does
#' not move. Ignored unless \code{algorithm = "discount"}, and a vector in
#' either produces one sub-model per value. See
#' \code{\link{create_vecxsubmodel}}.
#' @param algorithm algorithm that should be used for posterior simulation. If \code{NULL}
#' (default), standard algorithms will be used. See 'Details' for available
#' non-standard options.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' @param thin an integer thinning interval of the sampler (defaults to 1). After
#' the burn-in the sampler keeps the last of every \code{thin} draws, so it runs
#' \code{burnin + iterations * thin} draws and still keeps \code{iterations}.
#' See \code{\link{create_vecxsubmodel}}.
#' @param ... further arguments passed to or from other methods.
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
#'                         p_endogen = 1, p_exogen = 2,
#'                         global = "poil",
#'                         s = 2,
#'                         r = 1,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' 
#' @export
add_submodels.gvecmodel <- function(object,
                                    endogen = NULL,
                                    p_endogen = 1,
                                    exogen = NULL,
                                    p_exogen = 2,
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
                                    volatility_break = NULL,
                                    volatility_scales = 3,
                                    student_t = FALSE,
                                    delta_beta = 1,
                                    delta_sigma = 1,
                                    iterations = 10000,
                                    burnin = 2000,
                                    thin = 1,
                                    ...){
  # Student t errors are a VAR feature: BayesTS draws the scale of a period in
  # VarNormalWishart and in no error correction sampler. Asked for here the
  # argument would otherwise fall into '...' and be dropped, leaving a model
  # that is not the one the caller asked for and says nothing about it.
  if (!missing(student_t) && !identical(student_t, FALSE)) {
    stop("Student t errors are not available for an error correction ",
         "sub-model: BayesTS draws the scale of each period in the VAR with a ",
         "Wishart prior alone. Build the global model with ",
         "create_gvarmodel() to ask for them, or drop 'student_t'.",
         call. = FALSE)
  }
  
  submodels <- unique(object[["global"]][["index"]][, "submodel"])
  
  for (i in submodels) {
    object[["submodels"]][[i]] <- create_vecxsubmodel(object, submodel = i,
                                                      endogen = endogen, p_endogen = p_endogen,
                                                      exogen = exogen, p_exogen = p_exogen,
                                                      global = global, s = s,
                                                      r = r,
                                                      const = const, trend = trend, seasonal = seasonal,
                                                      structural = structural,
                                                      tvp = tvp, error = error, varsel = varsel,
                                                      algorithm = algorithm,
                                                      volatility_break = volatility_break,
                                                      volatility_scales = volatility_scales,
                                                      delta_beta = delta_beta,
                                                      delta_sigma = delta_sigma,
                                                      iterations = iterations, burnin = burnin,
                                                      thin = thin)
  }
  
  
    
  return(object)
}
