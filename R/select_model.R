#' Model Selection
#' 
#' Selects the best country-specific model of a GVAR model.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{estimate_gvar}}.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are "AIC", "BIC" (default) and "HQ".
#' @param select a character specifying how the best country model is selected.
#' If \code{"order"} (default), the country model with the overall minimum value of the
#' specified information criterion per country is selected. If \code{"rank"}, the
#' function selects the model, after which the selected information crition increases
#' for the first time.
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
#' # Add priors
#' object <- add_priors(object)
#' 
#' # Estimate GVAR model
#' gvar_est <- estimate_gvar(object, iterations = 100, burnin = 10, thin = 2)
#' # Note that the number of iterations and burnin should be much higher, e.g., 30000.
#' 
#' # Select models
#' gvar_est <- select_model(gvar_est)
#' 
#' @export
select_model <- function(object, ic = "BIC", select = "order") {
 
  criteria <- teststats(object, ic = ic, show_min = TRUE, select = select)
  
  pos_final <- which(criteria[, "Minimum"])
  result <- NULL
  n_final <- NULL
  for (i in pos_final){
    result <- c(result, list(object[[i]]))
    n_final <- c(n_final, names(object)[i])
  }
  names(result) <- n_final
  rm(n_final)
  
  class(result) <- append("bgvarest", class(result))
  return(result)
}