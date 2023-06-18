#' Specifications of the Submodels of a GVAR Model
#'
#' Obtains the model specification of the country-specific VARX or 
#' VECX models of a GVAR model.
#'
#' @param object an object of class \code{"bgvarest"} or \code{"bgvecest"}, usually,
#' a result of a call to \code{\link{draw_posterior.gvarsubmodels}} or 
#' \code{\link{draw_posterior.gvecsubmodels}}, respectively.
#'
#' @return A data frame.
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
#'                  domestic = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  global = list(variables = c("poil"), lags = 1),
#'                  deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                  iterations = 10,
#'                  burnin = 10)
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
#' # Obtain specifications
#' get_submodel_specifications(object)
#' 
#'
#' @export
get_submodel_specifications <- function(object){
  
  n_models <- length(object)
  
  result <- data.frame(ctry = rep(NA, n_models),
                       type = rep(NA, n_models),
                       r = rep(NA, n_models),
                       var_domestic = rep(NA, n_models),
                       lag_domestic = rep(NA, n_models),
                       var_foreign = rep(NA, n_models),
                       lag_foreign = rep(NA, n_models),
                       var_global = rep(NA, n_models),
                       lag_global = rep(NA, n_models),
                       ssvs = rep(NA, n_models),
                       bvs = rep(NA, n_models),
                       stringsAsFactors = FALSE)
  
  for (i in 1:n_models) {
    result[i, "ctry"] <- names(object)[i]
    type <- object[[i]][["model"]][["type"]]
    if (object[[i]][["model"]][["structural"]]) {
     type <- paste0("S", type) 
    }
    if (object[[i]][["model"]][["sv"]]) {
      type <- paste0("SV-", type) 
    }
    if (object[[i]][["model"]][["tvp"]]) {
      type <- paste0("TVP-", type) 
    }
    result[i, "type"] <- type
    rm(type)
    result[i, "var_domestic"] <- paste(object[[i]][["model"]][["domestic"]][["variables"]], collapse = ", ")
    result[i, "lag_domestic"] <- object[[i]][["model"]][["domestic"]][["lags"]]
    result[i, "var_foreign"] <- paste(object[[i]][["model"]][["foreign"]][["variables"]], collapse = ", ")
    result[i, "lag_foreign"] <- object[[i]][["model"]][["foreign"]][["lags"]]
    if (!is.null(object[[i]][["model"]][["global"]])) {
      result[i, "var_global"] <- paste(object[[i]][["model"]][["global"]][["variables"]], collapse = ", ")
      result[i, "lag_global"] <- object[[i]][["model"]][["global"]][["lags"]] 
    }
    if (!is.null(object[[i]][["model"]][["rank"]])) {
      result[i, "r"] <- object[[i]][["model"]][["rank"]] 
    }
    if (!is.null(object[[i]][["model"]][["varselect"]])) {
      result[i, "varselect"] <- object[[i]][["model"]][["varselect"]] 
    }
  }
  
  result <- result[, !unlist(lapply(result, function(x) {all(is.na(x))}))]
  
  return(result)
}
