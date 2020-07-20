#' Create Country Model Specifications
#' 
#' Produces a list of model specifications for each country in a GVAR model.
#' 
#' @param country_data a named list of time-series objects of country-specific data.
#' @param global_data a named time-series object of global data.
#' @param countries a character vector of country names, which should enter the GVAR model.
#' The names must correspond to the names of the elements in \code{country_data}. If
#' \code{NULL}, (default) all countries are considered.
#' @param variables a character vector specifying the considered variables. If \code{NULL},
#' all available variables will be used.
#' @param type a character specifying if the country models are estimated as vector
#' autoregressive "VAR" (default) or error correction "VEC" models.
#' @param p_domestic an integer or integer vector of the lag order of domestic variables in the VAR. See also 'Details'.
#' @param p_foreign an integer or integer vector of the lag order of foreign variables in the VAR. See also 'Details'.
#' @param s an integer or integer vector of the lag order of global variables in the VAR. See also 'Details'.
#' @param r an integer or vector for the cointegration rank. See also 'Details'.
#' 
#' @details
#' 
#' \code{p_domestic}, \code{p_foreign} and \code{s} refer to the lag order of the VAR model If \code{type = "VEC"},
#' the function \code{\link{create_models}} will automatically reduce each lag order by 1. If a vector is provided,
#' \code{\link{create_models}} will produce a distinct model for all possible specifications.
#' 
#' The argument \code{rank} is only used when \code{type = "VEC"}. Otherwise, the rank is set to zero.
#' If a vector is provided, \code{\link{create_models}} will produce a distinct model for all possible specifications.
#' 
#' @return The function produces a list of model specifications for each country,
#' where each element consists of the following elements:
#' \item{domestic}{a list of domestic variable names and their lag order.}
#' \item{foreign}{a list of foreign variable names and their lag order.}
#' \item{global}{a list of global variable names and their lag order (optional).}
#' \item{cointegration}{if \code{type = "VEC"}, a vector of rank specifications for
#' the cointegration matrix.}
#' \item{type}{a character specifying the type of the model, either \code{"VAR"} or \code{"VEC"}.}
#' 
#' @examples
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' global_data <- gvar2016$global_data
#' region_weights <- gvar2016$region_weights
#' weight_data <- gvar2016$weight_data
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp"),
#'                                      p_domestic = 2,
#'                                      p_foreign = 2,
#'                                      type = "VAR")
#' 
#' @export
create_specifications <- function(country_data, global_data = NULL,
                                  countries = NULL, variables = NULL,
                                  type = "VAR", p_domestic = 2, p_foreign = 2, s = 2, r = NULL){
  
  if (type == "VEC" & (0 %in% c(p_domestic, p_foreign, s))) {
    stop("The lag order for VEC models must be at least 1.")
  }
  
  if (!type %in% c("VAR", "VEC")) {"Argument 'type' must be either 'VAR' or 'VEC'."}
  
  #### Basic data checks ####
  # Check if country data is a list
  if (class(country_data) != "list") {
    stop("Country data must have class 'list'.")
  }
  # Check if data in country_data are of class ts
  if(sum(unlist(lapply(country_data, class)) == "ts") != length(country_data)) {
    stop("Elements in 'country_data' must be of class 'ts'.")
  }
  # Check if global data is of class ts
  if (!is.null(global_data)) {
    if(!"ts" %in% class(global_data)) {
      stop("'global_data' must be of class 'ts'.")
    }
  }
  
  # Check if data are named
  if(is.null(names(country_data))){
    stop("'country_data' must be a named list. Please provide a name for each country.")
  }
  if(!is.null(global_data) & is.null(dimnames(global_data)[[2]])){
    stop("'global_data' must be named. Please, provide a name for each global variable.")
  }
  
  #### Prepare specifications ####
  if (any(p_foreign < 1)) {
    stop("Lag of foreign star variables in VARX model must be at least 1.")
  }
  s.vars <- unique(unlist(lapply(country_data, function(x) {return(dimnames(x)[[2]])})))
  if (!is.null(variables)){
    s.vars <- s.vars[which(is.element(s.vars, variables))]
  }
  
  global <- FALSE
  if (is.null(global_data)){
    g.vars <- NA
  } else {
    g.vars <- dimnames(global_data)[[2]]
    if (!is.null(variables)) {
      g.used <- which(is.element(g.vars, variables)) 
    } else {
      g.used <- 1:length(g.vars)
    }
    if (length(g.used) > 0){
      g.vars <- g.vars[g.used]
      global <- TRUE
    }
  }
  
  specs <- list()
  specs.names <- NULL
  for (i in 1:length(country_data)){
    d.vars <- dimnames(country_data[[i]])[[2]]
    if (!is.null(variables)){
      d.vars <- d.vars[which(is.element(d.vars, variables))]
    }
    specs <- c(specs, list(list(domestic = list("variables" = d.vars, "lags" = p_domestic),
                                foreign = list("variables" = s.vars, "lags" = p_foreign))))
    if (global) {
      specs[[i]]$global = list("variables" = g.vars, "lags" = s)
    }
    if (type == "VEC") {
      if (is.null(r)) {
        specs[[i]]$cointegration = list(rank = 0:length(d.vars))
      } else {
        specs[[i]]$cointegration = list(rank = r)
      }
    }
    specs[[i]]$type <- type
    specs.names <- c(specs.names, names(country_data)[i])
  }
  names(specs) <- specs.names
  for (i in 1:length(specs)){
    specs[[i]]$domestic$variables <- specs[[i]]$domestic$variables[is.element(specs[[i]]$domestic$variables, dimnames(country_data[[i]])[[2]])]
  }
  if (!is.null(countries)){
    if (any(!is.element(countries, names(country_data)))) {
      stop(paste("There are no country data available for ", paste(countries[which(!is.element(countries, names(country_data)))], collapse = ", "), ".", sep = ""))
    }
    result <- NULL
    names.result <- NULL
    for (i in countries){
      result <- c(result, list(specs[[i]]))
      names.result <- c(names.result, i)
    }
    names(result) <- names.result
  } else {
    result <- specs
  }
  return(result)
}