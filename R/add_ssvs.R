#' Use Stochastic Search Variable Selection
#'
#' Adds specifictions for the use of stochastic Search Variable Selection as proposed
#' by George et al. (2008) to a list of country models, which was produced by the
#' function \code{\link{create_models}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_models}}.
#' @param exclude_deterministics logical. If \code{TRUE} (default), only non-deterministic terms
#' are subject to BVS during the estimation.
#' @param threshold a numeric between 0 and 1 specifiying the minimum posterior inclusion probability
#' of a variable to be included in the final model. If \code{NULL} (default), SSVS will run until the
#' end of the Gibbs sampler and no final model is chosen.
#' @param draws an integer specifying the amount of post-burn-in draws used to calculate posterior inclusion
#' probabilities.
#' 
#' @return A list of country models.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \url{https://doi.org/10.1016/j.jeconom.2007.08.017}
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
#' # Generate EA area region with 2 year, rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "MX", "GB"),
#'                                      p_domestic = 1,
#'                                      p_foreign = 1,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         model_specs = model_specs)
#' 
#' # Add SSVS
#' object <- add_ssvs(object)
#' 
#' @export
add_ssvs <- function(object, exclude_deterministics = TRUE,
                    threshold = NULL, draws = NULL){
  
  if (!is.null(threshold) & is.null(draws)) {
    stop("Argument 'draws' must be specified, if 'threshold' is specified.")
  }
  
  if (is.null(threshold) & !is.null(draws)) {
    stop("Argument 'threshold' must be specified, if 'draws' is specified.")
  }
  
  for (i in 1:length(object)) {
    if (!is.null(object[[i]]$variable_selection)) {
      stop("Variable selection method is already specified.")
    }
    
    object[[i]]$model$variable_selection$type <- "ssvs"
    object[[i]]$model$variable_selection$exclude_deterministics <- exclude_deterministics
    if (!is.null(threshold)) {
      object[[i]]$model$variable_selection$draws <- draws
      object[[i]]$model$variable_selection$threshold <- threshold
    }
  }
  return(object)
}