#' Contemporaneous Coefficients
#' 
#' Calculates statistics for coefficient draws of contemporaneous foreign variables
#' and their domestic counterparts for the country models of a Bayesian GVAR model.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{estimate_gvar}}.
#' @param type a character specifying the type of statistic that should be calculated.
#' Possible option are \code{"mean"} (default), \code{"median"} and \code{"sd"}.
#' @param t a numeric specifying the time index of the coefficients. Only applicable if the country models were
#' estimated with time-varying parameters or stochastic volatility.
#' If \code{NULL} the latest or only available period will be used.
#' @param all a logical. If \code{TRUE}, the contemporaneous effects of all estimated models will be printed. Otherwise, only
#' those models will be considered, which enter the global model.
#' @param variables a character vector specifying the variables, for which contemporaneous coefficients
#' should be extracted.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are \code{"AIC"}, \code{"BIC"} (default) and \code{"HQ"}.
#' @param select a character specifying how the best country model is selected.
#' If \code{"order"} (default), the country model with the overall minimum value of the
#' specified information criterion per country is selected. If \code{"rank"}, the
#' function selects the model, after which the selected information crition increases
#' for the first time.
#' 
#' @return Data frame.
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
#' # Generate EA area region with 3 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 3,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 3 year rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 3,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "GB"),
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
#' gvar_est <- estimate_gvar(object, iterations = 50, burnin = 10, thin = 4)
#' # Note that the number of iterations and burnin should be much higher.
#' 
#' contemp_coeffs(gvar_est)
#' 
#' @export
contemp_coeffs <- function(object, type = "mean", t = NULL, all = FALSE,
                                 variables = NULL, ic = "BIC", select = "order"){
  
  data <- select_model(object, ic = ic, select = select)
  
  # Get available variables
  countries <- names(data)
  vars_domestic <- unique(unlist(lapply(data, function(x){return(x$model$domestic$variables)})))
  vars_foreign <- unique(unlist(lapply(data, function(x){return(x$model$foreign$variables)})))
  if (is.null(variables)) {
    variables <- vars_foreign[which(is.element(vars_foreign, vars_domestic))]
    variables <- variables[order(variables)] 
  }
  
  result <- as.data.frame(matrix(NA, length(countries), length(variables) + 1))
  names(result) <- c("Country", variables)
  result[, "Country"] <- countries
  
  tt <- t
  
  for (i in countries){
    pos_i <- which(is.element(countries, i)) # Position in result table
    var_i_domestic <- data[[i]]$model$domestic$variables # Domestic variables
    n_i_domestic <- length(var_i_domestic)
    var_i_foreign <- data[[i]]$model$foreign$variables # Foreign star variables
    n_i_foreign <- length(var_i_foreign)
    use <- var_i_domestic[which(is.element(var_i_domestic, var_i_foreign))] # Find valid cases
    tvp <- FALSE
    if (length(use) > 0) {
      # Check if tvp
      if (class(data[[i]]$coefficients$a_foreign) == "list") {
        tvp <- TRUE
        if (is.null(t)) {
          tt <- length(data[[i]]$coefficients$a_foreign)
        }
      }
      
      # Extract contemporaneous effects
      if (tvp) {
        contemp <- data[[i]]$coefficients$a_foreign[[tt]][, 1:(n_i_domestic * n_i_foreign)]
      } else {
        contemp <- data[[i]]$coefficients$a_foreign[, 1:(n_i_domestic * n_i_foreign)] 
      }
      
      if (type == "mean") {
        contemp <- matrix(apply(contemp, 2, mean), n_i_domestic)
      }
      if (type == "median") {
        contemp <- matrix(apply(contemp, 2, stats::median), n_i_domestic)
      }
      if (type == "sd") {
        contemp <- matrix(apply(contemp, 2, stats::sd), n_i_domestic)
      }
      dimnames(contemp) <- list(var_i_domestic, var_i_foreign)
      for (j in use) {
        result[pos_i, which(names(result) == j)] <- contemp[j, j]
      }
    }
  }
  
  # Omit redundant columns
  na_cols <- which(apply(result, 2, function(x) {all(is.na(x))}))
  if (length(na_cols) > 0) {
    result <- result[, -na_cols]
  }
  rm(na_cols)
  row.names(result) <- NULL # Reset row names
  return(result)
}