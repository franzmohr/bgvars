#' Submodel Selection
#' 
#' Selects the best country-specific submodel of a GVAR model.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{draw_posterior}}.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are "AIC", "BIC" (default) and "HQ".
#' @param select a character specifying how the best country model is selected.
#' If \code{"order"} (default), the country model with the overall minimum value of the
#' specified information criterion per country is selected. If \code{"rank"}, the
#' function selects the model, after which the selected information criterion increases
#' for the first time.
#' @param teststats optional. An object of class "ctryvartest", usually,
#' the result of a call to \code{\link{submodel_test_statistics}}.
#' 
#' @examples
#' # Load data
#' data("gvar2019")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2019$country_data,
#'              weight_data = gvar2019$weight_data,
#'              region_weights = gvar2019$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2019$global_data
#' 
#' # Difference series to make them stationary
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
#' global_data <- diff_variables(global_data, multi = 100)
#' 
#' # Create time varying weights
#' weight_data <- create_weights(weight_data, period = 3, country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(
#'                  country_data = country_data,
#'                  global_data = global_data,
#'                  countries = c("US", "JP", "CA", "NO", "GB", "EA"), 
#'                  domestic = list(variables = c("y", "Dp", "r"), lags = 1:2),
#'                  foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  global = list(variables = c("poil"), lags = 1),
#'                  deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                  iterations = 4,
#'                  burnin = 2)
#' # Note that the number of iterations and burnin draws should be much higher!
#'                                      
#' # Overwrite country-specific specifications
#' model_specs[["US"]][["domestic"]][["variables"]] <- c("y", "Dp", "r")
#' model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")
#' 
#' # Create estimation objects
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#' 
#' # Add priors
#' models_with_priors <- add_priors(country_models,
#'                                  coef = list(v_i = 1 / 9, v_i_det = 1 / 10),
#'                                  sigma = list(df = 3, scale = .0001))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(models_with_priors)
#' 
#' # Select final submodels
#' object <- select_submodels(object, ic = "BIC", select = "order")
#' 
#' 
#' @export
select_submodels <- function(object, ic = "BIC", select = "order", teststats = NULL) {
  
  obj_class <- class(object)
  
  if (!ic %in% c("AIC", "BIC", "HQ")) {
    stop("Invalid specification of argument 'ic'.")
  }
  
  if (length(select) > 1) {
    stop("Invalid specification of argument 'select'.")
  }
  
  if (!select %in% c("order", "rank")) {
    stop("Invalid specification of argument 'select'.")
  }
  
  # Obtain test statistics
  if (is.null(teststats)) {
    criteria <- submodel_test_statistics(object)
  } else {
    # Basic checks
    if (!"teststats.bgvarest" %in% class(teststats)) {
      stop("If provided, argument 'teststats' must be of class teststats.bgvarest")
    }
    criteria <- teststats
  }
  
  names_object <- names(object)
  
  crit <- list()
  ctry <- unique(names(object))
  for (i in 1:length(ctry)) {
    crit[[i]] <- criteria[[i]][["teststats"]]
  }
  names(crit) <- ctry
  criteria <- crit
  
  # Select final country models
  result <- NULL
  for (i in 1:length(criteria)) {
    # Order selection
    if (select == "order") {
      sub_pos <- which.min(criteria[[i]][, ic])
    }
    # Add rank selection here
    if (select == "rank") {
      crit <- criteria[[i]][, ic]
      sub_pos <- 1
      if (length(crit) > 1) {
        while (crit[sub_pos + 1] < crit[sub_pos]) {
          sub_pos <- sub_pos + 1
        }
      }
    }
    pos <- which(names_object == names(criteria)[i])[sub_pos]
    result[[i]] <- object[[pos]]
  }
  names(result) <- names(criteria)
  
  class(result) <- obj_class
  return(result)
}