

.input_check_countries <- function(countries, country_data) {
  if (!is.null(countries)) {
    if (length(countries) < 2) {
      stop("Please, specifiy at least two countries.")
    }
    if (!all(countries %in% names(country_data))) {
      stop("Not all countries specified in argument 'countries' are available in argument 'country_data'.")
    }
  }
}

.input_check_country_data <- function(country_data) {
  if (!"list" %in% class(country_data)) {
    stop("Argument 'country_data' must be of class 'list'.")
  }
  # Check if data in country_data are of class ts
  if(any(!unlist(lapply(country_data, function(x) {"ts" %in% class(x)})))) {
    stop("Elements in 'country_data' must be of class 'ts'.")
  }
  if(is.null(names(country_data))){
    stop("'country_data' must be a named list.")
  } 
}

.input_check_global_data <- function(global_data) {
  if (!is.null(global_data)) {
    if(!"ts" %in% class(global_data)) {
      stop("'global_data' must be of class 'ts'.")
    }
    
    if(is.null(dimnames(global_data)[[2]])){
      stop("'global_data' must be named. Please, provide a name for each global variable.")
    }
  }
}

.input_check_model_specs <- function(x, name) {
  if (!"list" %in% class(x)) {
    stop(paste0("Argument '", name, "' must be a named list."))
  }
  if (is.null(x[["lags"]])) {
    stop(paste0("Element 'lags' is missing in argument '", name, "'."))
  }
  if (any(x$lags < 0)) {
    stop(paste0("Values of element 'lags' in argument '", name, "' must be at least 0."))
  }
}

.input_check_structural <- function(structural, countries) {
  if (!"character" %in% class(structural)) {
    stop("Argument 'structural' must be of class character.")
  }
  
  if (!is.null(countries)) {
    if (!all(structural %in% countries)) {
      stop("Not all countries in argument 'structural' are available.")
    }
  }
}

.input_check_weight_data <- function(weight_data) {
  # Check class
  if (!any(c("array", "matrix") %in% class(weight_data))) {
    stop("Argument 'weight_data' must be an array or a matrix.")
  }
  # Check if weights are named
  if (is.null(dimnames(weight_data)[[1]]) | is.null(dimnames(weight_data)[[2]])) {
    stop("Rows and columns of 'weight_data' must both be named.")
  }
}