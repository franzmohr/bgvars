#' Country Model Test Statistics
#' 
#' Extracts the specifications of the estimated country-specific models of a GVAR model.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{estimate_gvar}}.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are "AIC", "BIC" (default) and "HQ".
#' @param show_min a logical indicating whether a column should be added, which indicates
#' the best country model for the specified information criterion.
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
#' # Generate EA area region with 2 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year rolling window averages
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
#' # Get model specifications
#' gvar_est <- teststats(gvar_est)
#' 
#' @export
teststats <- function(object, ic = "BIC", show_min = TRUE, select = "order") {
  
  if (!ic %in% c("AIC", "BIC", "HQ")) {
    stop("Information criterion not available. Only 'AIC', 'BIC' or 'HQ' are available.")
  }
  
  if (!select %in% c("rank", "order")) {
    stop("Argument 'select' must be either 'order' or 'rank'.")
  }
  
  criteria <- data.frame("Country" = names(object), "Rank" = NA,
                         "p_domestic" = NA, "p_foreign" = NA, "s" = NA,
                         "AIC" = NA, "BIC" = NA, "HQ" = NA, stringsAsFactors = FALSE)

  for (i in 1:length(object)) {
    if (!is.null(object[[i]]$coefficients)) {
      criteria[i, c("AIC", "BIC", "HQ")] <- object[[i]]$teststats[c("AIC", "BIC", "HQ")]
      if (!is.null(object[[i]]$model$cointegration$rank)) {
        criteria[i, "Rank"] <- object[[i]]$model$cointegration$rank
      }
      criteria[i, "p_domestic"] <- object[[i]]$model$domestic$lags
      criteria[i, "p_foreign"] <- object[[i]]$model$foreign$lags
      if (!is.null(object[[i]]$model$global)) {
        criteria[i, "s"] <- object[[i]]$model$global$lags 
      } 
    }
  }
  
  if (show_min) {
    criteria <- cbind(criteria, "Minimum" = FALSE)
    for (i in unique(names(object))) {
      c_pos <- which(criteria[, "Country"] == i)
      if (length(c_pos) > 1) {
        c_temp <- criteria[c_pos, ic]
        if (select == "rank") {
          k <- 1
          while (!is.na(c_temp[k + 1]) & c_temp[k] > c_temp[k + 1] & k < length(c_temp)) {
            k <- k + 1
          }
        }
        if (select == "order") {
          k <- which.min(c_temp)
        }
        criteria[c_pos[k], "Minimum"] <- TRUE 
      } else {
        criteria[c_pos, "Minimum"] <- TRUE 
      }
    }
  }
  
  na_cols <- which(apply(criteria, 2, function(x) {all(is.na(x))}))
  if (length(na_cols) > 0) {
    criteria <- criteria[, -na_cols]
  }
  rm(na_cols)
  row.names(criteria) <- NULL
  
  return(criteria)
}