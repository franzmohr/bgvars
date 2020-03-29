#' Add Deterministic Terms to Country Models
#'
#' Adds deterministic terms to a list of country models, which was produced by
#' the function \code{\link{create_models}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_models}}.
#' @param const either logical for VAR models or a character specifying whether a constant term
#' enters the error correction term (\code{"restricted"}) or the non-cointegration term (\code{"unrestricted"})
#' of a VEC model. If `NULL` (default) no constant term will be added.
#' @param trend either logical for VAR models or a character specifying whether a trend term
#' enters the error correction term (\code{"restricted"}) or the non-cointegration term (\code{"unrestricted"})
#' of a VEC model. If `NULL` (default) no trend term will be added.
#' @param seasonal either logical for VAR models or a character specifying whether seasonal dummies
#' enter the error correction term (\code{"restricted"}) or the non-cointegreation term (\code{"unrestricted"}).
#' If `NULL` (default) no seasonal terms will be added. The amount of dummy variables depends
#' on the frequency of the time series data in each list element.
#'  
#' @return A list of country models.
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
#' # Add deterministics
#' object <- add_deterministics(object, const = TRUE)
#' 
#' # Add priors
#' object <- add_priors(object)
#' 
#' 
#' @export
add_deterministics <- function(object, const = NULL, trend = NULL, seasonal = NULL){
  
  # Check if priors have already been generated
  if (any(unlist(lapply(object, function(x){!is.null(x$priors)})))) {
    stop("Priors have already been generated. Please use 'add_deterministics' before 'add_priors'.")
  }
  
  for (i in 1:length(object)) {
    
    tt <- nrow(object[[i]]$data$x_domestic)
    temp_tsp <- stats::tsp(object[[i]]$data$x_domestic)
    
    if (object[[i]]$model$type == "VAR") {
      
      if (is.null(const)) {
        const <- FALSE
      } else {
        if (class(const) != "logical") {
          stop("Argument 'const' must be of class 'logical' for VAR models.")
        }
      }
      
      if (is.null(trend)) {
        trend <- FALSE
      } else {
        if (class(trend) != "logical") {
          stop("Argument 'trend' must be of class 'logical' for VAR models.")
        }
      }

      if (is.null(seasonal)) {
        seasonal <- FALSE
      } else {
        if (class(seasonal) != "logical") {
          stop("Argument 'seasonal' must be of class 'logical' for VAR models.")
        }
      }
      
      temp <- NULL
      temp_names <- NULL
      if (const) {
        temp <- cbind(temp, rep(1, tt))
        temp_names <- append(temp_names, "const")
      }
      if (trend) {
        temp <- cbind(temp, 1:tt)
        temp_names <- append(temp_names, "trend")
      }
      if (seasonal) {
        freq <- stats::frequency(object[[i]]$data$x_domestic)
        if (freq == 1) {
          warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
        } else {
          pos <- which(floor(stats::time(object[[i]]$data$x_domestic)) == stats::time(object[[i]]$data$x_domestic))[1]
          pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
          for (j in 1:(freq - 1)) {
            s_temp <- rep(0, freq)
            s_temp[pos[j]] <- 1
            temp <- cbind(temp, rep(s_temp, length.out = tt))
            temp_names <- c(temp_names, paste("season.", j, sep = ""))
          }
        }
      }
      
      temp <- stats::ts(as.matrix(temp), class = c("mts", "ts", "matrix"))
      stats::tsp(temp) <- temp_tsp
      dimnames(temp)[[2]] <- temp_names
    } else {
      if (!is.null(const) & class(const) == "logical") {
        stop("Argument 'const' must not be logical for VEC models.")
      }
      if (!is.null(trend) & class(trend) == "logical") {
        stop("Argument 'trend' must not be logical for VEC models.")
      }
      if (!is.null(seasonal) & class(seasonal) == "logical") {
        stop("Argument 'seasonal' must not be logical for VEC models.")
      }
      temp_res <- NULL
      temp_unres <- NULL
      temp_names_res <- NULL
      temp_names_unres <- NULL
      
      if (!is.null(const)) {
        if (const == "restricted") {
          temp_res <- cbind(temp_res, rep(1, tt))
          temp_names_res <- append(temp_names_res, "const")
        }
        if (const == "unrestricted") {
          temp_unres <- cbind(temp_unres, rep(1, tt))
          temp_names_unres <- append(temp_names_unres, "const")
        } 
      }
      
      if (!is.null(trend)) {
        if (trend == "restricted") {
          temp_res <- cbind(temp_res, 1:tt)
          temp_names_res <- append(temp_names_res, "trend")
        }
        if (trend == "unrestricted") {
          temp_unres <- cbind(temp_unres, 1:tt)
          temp_names_unres <- append(temp_names_unres, "trend")
        }
      }
      
      if (!is.null(seasonal)) {
        if (seasonal %in% c("restricted", "unrestricted")) {
          seas <- NULL
          freq <- stats::frequency(object[[i]]$data$x_domestic)
          if (freq == 1) {
            warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
          } else {
            pos <- which(floor(stats::time(object[[i]]$data$x_domestic)) == stats::time(object[[i]]$data$x_domestic))[1]
            pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
            for (j in 1:(freq - 1)) {
              s_temp <- rep(0, freq)
              s_temp[pos[j]] <- 1
              if (seasonal == "restricted") {
                temp_res <- cbind(temp_res, rep(s_temp, length.out = tt))
                temp_names_res <- c(temp_names_res, paste("season.", j, sep = "")) 
              }
              if (seasonal == "unrestricted") {
                temp_unres <- cbind(temp_unres, rep(s_temp, length.out = tt))
                temp_names_unres <- c(temp_names_unres, paste("season.", j, sep = ""))
              }
            }
          }
        }
      }
      
      temp <- list()
      if (!is.null(temp_res)) {
        temp_res <- stats::ts(as.matrix(temp_res), class = c("mts", "ts", "matrix"))
        stats::tsp(temp_res) <- temp_tsp
        dimnames(temp_res)[[2]] <- temp_names_res
        temp$restricted <- temp_res
        rm(temp_res)
      }
      if (!is.null(temp_unres)) {
        temp_unres <- stats::ts(as.matrix(temp_unres), class = c("mts", "ts", "matrix"))
        stats::tsp(temp_unres) <- temp_tsp
        dimnames(temp_unres)[[2]] <- temp_names_unres
        temp$unrestricted <- temp_unres
        rm(temp_unres)
      }
    }
    object[[i]]$data$deterministics <- temp
    rm(list = c("temp", "temp_tsp"))
  }
  return(object)
}