#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class "bvarprd".
#' 
#' @param x an object of class "bvarprd", usually, a result of a call to \code{\link{predict.bgvar}}.
#' @param variable a character vector containing a country name as its first and a variable name
#' as its second element. If \code{NULL} (default) all series are plotted.
#' @param ... further graphical parameters.
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
#' # Note that the number of iterations and burnin should be much higher.
#' 
#' # Solve GVAR
#' gvar_solved <- solve_gvar(gvar_est)
#' 
#' # Forecasts
#' gvar_predict <- predict(gvar_solved)
#'
#' # Plot Forecasts
#' plot(gvar_predict, variable = c("US", "Dp"))
#' 
#' @export
plot.bgvarprd <- function(x, variable = NULL, ...) {
  
  tt <- nrow(x[["y"]])
  if (is.null(variable)) {
    pos <- 1:ncol(x[["y"]])
  } else {
    pos <- which(x[["index"]][, "country"] == variable[1] & x$index[, "variable"] == variable[2])
  }
  for (i in pos) {
    n_ahead <- nrow(x[["fcst"]][[i]])
    temp <- matrix(NA, tt + n_ahead, 4)
    temp[1:tt, 1] <- x[["y"]][, i]
    temp[tt, 2:4] <- x[["y"]][tt, i]
    temp[(tt + 1):(tt + n_ahead), 2:4] <- x[["fcst"]][[i]]
    tsp_temp <- stats::tsp(x[["fcst"]][[i]])
    temp <- stats::ts(temp, end = tsp_temp[2], frequency = tsp_temp[3])
    rm(tsp_temp)
    title_temp <- paste(x[["index"]][i, "country"], ": ", x[["index"]][i, "variable"], sep = "")
    stats::plot.ts(temp, plot.type = "single", lty = c(1, 2, 1, 2), main = title_temp, ylab = "")
  }
  
}