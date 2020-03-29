#' Use Bayesian Variable Selection
#'
#' Adds specifictions for the use of Bayesian variable selection (BVS) as proposed in Korobilis (2013) to
#' a list of country models, which was produced by the function \code{\link{create_models}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_models}}.
#' @param exclude_deterministics logical. If \code{TRUE} (default), only non-deterministic terms
#' are subject to BVS during the estimation.
#' @param threshold a numeric between 0 and 1 specifiying the minimum posterior inclusion probability
#' of a variable to be included in the final model. If \code{NULL} (default), BVS will run until the
#' end of the Gibbs sampler and no final model is chosen.
#' @param draws an integer specifying the amount of post-burn-in draws used to calculate posterior inclusion
#' probabilities.
#' 
#' @return A list of country models.
#' 
#' @references
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230.
#' 
#' @examples
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' global_data <- gvar2016$global_data
#' region_weights <- gvar2016$region_weights
#' weight_data <- gvar2016$weight_data
#' 
#' # Take first difference of all variables y and Dp
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"))
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("US", "JP", "CA", "MX", "GB"),
#'                                      p_domestic = 1,
#'                                      p_foreign = 1,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         model_specs = model_specs)
#' 
#' # Add BVS
#' object <- add_bvs(object)
#' 
#' @export
add_bvs <- function(object, exclude_deterministics = TRUE,
                    threshold = NULL, draws = NULL){
  
  if (!is.null(threshold) & is.null(draws)) {
    stop("Argument 'draws' must be specified, if 'threshold' is specified.")
  }
  
  if (is.null(threshold) & !is.null(draws)) {
    stop("Argument 'threshold' must be specified, if 'draws' is specified.")
  }
  
  for (i in 1:length(object)) {
    if (object[[i]]$model$variable_selection != "none") {
      stop("Variable selection method is already specified.")
    }
    
    object[[i]]$model$variable_selection$type <- "bvs"
    object[[i]]$model$variable_selection$exclude_deterministics <- exclude_deterministics
    if (!is.null(threshold)) {
      object[[i]]$model$variable_selection$draws <- draws
      object[[i]]$model$variable_selection$threshold <- threshold
    }
  }
  return(object)
}