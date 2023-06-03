#' Create Country Models
#'
#' Produces a list of country models for each model specification that should be estimated.
#'
#' @param country_data a named list of time-series objects of country-specific data.
#' @param weight_data a matrix or an array of weights, usually, the output of a call to \code{\link{create_weights}}.
#' @param global_data a named time-series object of global data.
#' @param model_specs a list of model specifications, usuall, the output of a call to \code{\link{create_specifications}}.
#' 
#' @return A list of country data and model specifications.
#' 
#' @examples
#' 
#' # Load gvar2016 dataset
#' data("gvar2016")
#'
#' # Data objects
#' country_data <- gvar2016[["country_data"]]
#' global_data <- gvar2016[["global_data"]]
#' region_weights <- gvar2016[["region_weights"]]
#' weight_data <- gvar2016[["weight_data"]]
#' 
#' # Obtain weights
#' weight_data <- create_weights(weight_data = weight_data, period = 3,
#'                               country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      countries = c("US", "JP", "CA", "NO", "GB"), 
#'                                      domestic = list(variables = c("y", "Dp", "r"), lags = 1:2),
#'                                      foreign = list(variables = c("y", "Dp", "r"), lags = 1:2),
#'                                      global = list(variables = c("poil"), lags = 1:2),
#'                                      deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                                      iterations = 10,
#'                                      burnin = 10,
#'                                      thin = 1)
#' 
#' # Create list element for each model
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#' 
#' @export
create_models <- function(country_data, weight_data, global_data = NULL,
                          model_specs = NULL){
  
  # rm(list = ls()[-which(ls() %in% c("country_data", "weight_data", "global_data", "model_specs"))]);
  
  # Check if weights are named
  if (is.null(dimnames(weight_data)[[1]]) | is.null(dimnames(weight_data)[[2]])) {
    stop("Rows and columns of 'weight_data' must both be named.")
  }
  
  #### Prepare weight matrices ####
  
  # Build skeletons for weight matrices
  countries <- names(model_specs)
  const_weight <- "matrix" %in% class(weight_data)
  if (const_weight) {
    temp_w <- matrix(0, length(countries), length(countries))
    dimnames(temp_w) <- list(countries, countries)
  } else {
    temp_w <- array(0, dim = c(length(countries), length(countries), dim(weight_data)[[3]]))
    dimnames(temp_w) <- list(countries, countries, dimnames(weight_data)[[3]])
  }
  
  # Recalculate weights for the used countries
  temp_c <- NULL
  for (i in countries){
    temp_c <- c(temp_c, list(country_data[[i]]))
    if (const_weight){
      temp_w[i,] <- weight_data[i, countries] / (1 - sum(weight_data[i, !is.element(dimnames(weight_data)[[2]], countries)]))
    } else {
      for (j in 1:dim(weight_data)[3]){
        temp_w[i,,j] <- weight_data[i, countries, j] / (1 - sum(weight_data[i, !is.element(dimnames(weight_data)[[2]], countries), j])) 
      }
    }
  }
  country_data <- temp_c; rm(temp_c)
  names(country_data) <- countries
  weight_data <- temp_w; rm(temp_w)
  
  # Check if all country time series have the same length ----
  n_c <- length(countries)
  t_max <- c()
  for (i in 1:n_c){
    if (any(apply(country_data[[i]], 2, function(x){any(is.na(x))}))){
      vars <- dimnames(country_data[[i]])[[2]][which(apply(country_data[[i]], 2, function(x){any(is.na(x))}))]
      country <- names(country_data)[i]
      stop(paste("In country ", country, " the variable(s) ", paste(vars, sep = ", "), " contain NA values.", sep = ""))
      rm(vars)
    }
    t_max <- append(t_max, dim(stats::na.omit(country_data[[i]]))[1])
  }
  if (min(t_max) != max(t_max)) {
    stop("Numbers of observations differ across countries. For now, this package requires the same number for each country.")
  } else {
    t <- t_max[1]
    rm(t_max)
  }
  
  # Check if weight matrix has the same order of names in the rows and columns
  if (any(dimnames(weight_data)[[1]] != dimnames(weight_data)[[2]])){
    stop("The order of country names in the rows and columns of the weight matrix is not the same.")
  }
  
  # Check if weight matrix contains data on all required countries
  if (any(!is.element(names(country_data), dimnames(weight_data)[[1]]))){
    stop("The weight matrix does not contain data for at least one country in the country sample or is named differently.")
  }
  
  # If weight matrix is time varying, check if the number of periods
  # is the same as in the country series
  if (!is.na(dim(weight_data)[3])){
    if (dim(weight_data)[3] > t) {
      warning("The weight array does not contain as much periods as the country sample. Trying to correct this.")
      dim_w <- as.numeric(dimnames(weight_data)[[3]])
      w_t <- stats::ts(dim_w, start = dim_w[1], frequency = stats::frequency(country_data[[1]]))
      weight_data <- weight_data[,,-which(!w_t %in% stats::time(country_data[[1]]))]
    }
  }
  
  # Reduce object with global data to relevant series  ----
  global_vars <- any(unlist(lapply(model_specs, function(x) {return(!is.null(x$global))})))
  if (global_vars){
    if (is.null(global_data)) {
      stop("The model specifications contain global variables, but argument 'global_data' was not specified.")
    }
    global_vars_temp <- unique(unlist(lapply(model_specs, function(x) {return(x$global$variables)})))
    tsp_temp <- stats::tsp(global_data)
    global_data <- stats::as.ts(as.matrix(global_data[, global_vars_temp]))
    dimnames(global_data)[[2]] <- global_vars_temp
    stats::tsp(global_data) <- tsp_temp
    
    # In case global data has different length than country data
    if (t < nrow(global_data)) {
      n_global <- NCOL(global_data)
      global_data <- na.omit(cbind(global_data, country_data[[1]]))[, 1:n_global]
      tsp_temp <- stats::tsp(global_data)
      global_data <- stats::as.ts(as.matrix(global_data))
      dimnames(global_data)[[2]] <- global_vars_temp
      stats::tsp(global_data) <- tsp_temp
    }
    
    if (t > nrow(global_data)) {
      pos_global <- 1:NCOL(global_data)
      for (i in names(country_data)) {
        temp_names <- dimnames(country_data[[i]])
        c_data <- na.omit(cbind(country_data[[i]], global_data))[, 1:length(temp_names[[2]])]
        tsp_temp <- stats::tsp(c_data)
        c_data <- stats::as.ts(as.matrix(c_data))
        dimnames(c_data)[[2]] <- temp_names[[2]]
        stats::tsp(c_data) <- tsp_temp
        country_data[[i]] <- c_data
        rm(c_data)
      }
      t <- nrow(country_data[[1]])
    }
    rm(tsp_temp)
  }
  
  if (!is.null(global_data) & !global_vars) {
    stop("Global data were provided in 'global_data', but they are not contained in the model specifications. Please adapt specifications, if you want to use global variables.")
  }
  rm(global_vars)
  
  # Generate final variable index ----
  # Get list of all variables
  avail.dom.vars <- unique(unlist(lapply(country_data, function(x) {return(dimnames(x)[[2]])}))) # Domestic variables
  if (!is.null(global_data)){
    avail.global.vars <- dimnames(global_data)[[2]] # Global variables
  }
  
  # Create index
  index <- data.frame() # Index
  X <- c() # Final data set of domestic variables (incl. global endogenous)
  X_names <- c()
  for (i in countries){
    temp <- c()
    temp_names <- c()
    for (j in model_specs[[i]]$domestic$variables){
      if (is.element(j, dimnames(country_data[[i]])[[2]])){
        temp <- cbind(temp, country_data[[i]][, j])
      } else {
        temp <- cbind(temp, global_data[, j])
      }
      temp_names <- append(temp_names, j)
    }
    X <- cbind(X, temp)
    n_temp <- length(temp_names)
    X_names <- c(X_names, temp_names)
    index <- rbind(index, data.frame("Country" = rep(i, n_temp),
                                     "Variable" = temp_names))
  }
  X <- stats::na.omit(X)
  dimnames(X)[[2]] <- X_names
  index[, 1] <- as.character(index[, 1])
  index[, 2] <- as.character(index[, 2])
  
  # Update foreign and global variables if global variable is endogenous in one country ----
  variables <- unique(index[, 2]) # All endogenous variables
  endo.g.var <- NULL
  for (i in countries){
    # Check if there is an endogenous global variable in a country specification
    if (any(is.element(model_specs[[i]]$global$variables, variables))){
      foreign <- model_specs[[i]]$foreign$variables
      global <- model_specs[[i]]$global$variables
      # Which global variable is endogenous
      pos.endo.g <- which(is.element(global, variables))
      # Endo global vars
      endo.global <- global[pos.endo.g]
      # Add to vector of endo global vars
      endo.g.var <- c(endo.g.var, endo.global)
      # If a global variable is endogenous in another country model
      # add it to the foreign variable vector
      if (!any(endo.global %in% model_specs[[i]]$domestic$variables)) {
        foreign <- c(foreign, endo.global)
        foreign <- unique(foreign)
      }
      # Remove from global variables
      global <- global[-pos.endo.g]
      if (length(global) == 0){
        model_specs[[i]]$global <- NULL
      } else {
        model_specs[[i]]$global$variables <- global
      }
      model_specs[[i]]$foreign$variables <- foreign
    }
  }
  
  # Update global data
  endo.g.var <- unique(endo.g.var)
  if (length(endo.g.var) > 0){
    if (length(endo.g.var) == length(dimnames(global_data)[[2]])){
      global_data <- NULL
    } else {
      global_data <- global_data[, -which(is.element(endo.g.var, dimnames(global_data)[[2]]))]
    }
  }
  
  # Trim endo and global data to same length ----
  global <- !is.null(global_data) # If a global variable is used
  if (global){
    k <- dim(X)[2]
    names_X <- dimnames(X)[[2]]
    names_g <- dimnames(global_data)[[2]]
    # Join endogenous and global data and omit NAs
    temp <- stats::na.omit(cbind(X, global_data))
    # Extract final endogenous variables
    X <- temp[, 1:k]
    # Extract final global variables
    X_global <- stats::as.ts(as.matrix(temp[, -(1:k)]))
    stats::tsp(X_global) <- stats::tsp(X)
    dimnames(X)[[2]] <- names_X
    dimnames(X_global)[[2]] <- names_g
    class(X_global) <- c("mts", "ts", "matrix")
    rm(list = c("temp", "names_X", "names_g"))
  } else {
    X_global <- NULL
  }
  rm(global_data)
  
  # Max length of global VAR model ---
  # Used so that each country model contains the same number of
  # observations so that information criteria are comparable
  max_lag <- c()
  for (i in countries) {
    max_lag <- append(max_lag, c(model_specs[[i]]$domestic$lags,
                                 model_specs[[i]]$foreign$lags))
  }
  if (global) {
    for (i in countries) {
      max_lag <- append(max_lag, model_specs[[i]]$global$lags)
    }
  }
  max_lag <- max(max_lag)
  for (i in countries){
    model_specs[[i]]$max_lag <- max_lag
  }
  
  # WEIGHT MATRICES ----
  
  # Generate weight matrices for each country ----
  tv <- !is.na(dim(weight_data)[3])
  
  # Check data availability for each country
  n_v <- length(variables)
  var.exists <- as.data.frame(matrix(NA, n_c, n_v))
  row.names(var.exists) <- countries
  names(var.exists) <- variables
  for (i in countries){
    for (j in variables){
      var.exists[i, j] <- is.element(j, index[index[,1] == i, 2])
    }
  }
  k <- dim(index)[1]
  
  W <- c()
  W.names <- c()
  for (i in countries){
    n.endo <- length(index[index[, 1] == i, 1])
    star.vars <- model_specs[[i]]$foreign$variables
    n.star <- length(star.vars)
    if (tv){
      W.i <- array(0, dim = c(n.endo + n.star, k, dim(weight_data)[3]))
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2], dimnames(weight_data)[[3]])
      for (l in 1:dim(weight_data)[3]) {
        W.i[1:n.endo, which(index[, 1] == i), l] <- diag(1, n.endo)
        for (j in star.vars){
          w.temp <- weight_data[i,,l]
          # Set all values to zero, where variable is not available
          w.temp[!var.exists[, j]] <- 0
          # Extract weights for available values
          exist.weights <- weight_data[i, var.exists[, j], l]
          # Make sum over weights of not available values
          non.exist.weights <- sum(weight_data[i, !var.exists[, j], l])
          # Recalculate weights
          w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
          # Add weights to country's weight matrix
          pos.x <- n.endo + which(star.vars == j)
          pos.y <- which(index[, 2] == j)
          W.i[pos.x, pos.y, l] <- w.temp[var.exists[, j]]
        } 
      }
    } else {
      W.i <- matrix(0, n.endo + n.star, k)
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2])
      W.i[1:n.endo, which(index[, 1] == i)] <- diag(1, n.endo)
      for (j in star.vars){
        w.temp <- weight_data[i,]
        # Set all values to zero, where variable is not available
        w.temp[!var.exists[, j]] <- 0
        # Extract weights for available values
        exist.weights <- weight_data[i, var.exists[, j]]
        # Make sum over weights of not available values
        non.exist.weights <- sum(weight_data[i, !var.exists[, j]])
        # Recalculate weights
        w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
        # Add weights to country's weight matrix
        pos.x <- n.endo + which(star.vars == j)
        pos.y <- which(index[, 2] == j)
        W.i[pos.x, pos.y] <- w.temp[var.exists[, j]]
      }
    }
    
    W <- c(W, list(W.i))
    W.names <- c(W.names, i)
  }
  names(W) <- W.names
  
  # Generate vector z = (x, x.star)' for each country ----
  X_temp <- t(X)
  Z <- NULL
  for (j in countries) {
    if ("matrix" %in% class(W[[j]])){
      x <- W[[j]] %*% X_temp
    } else {
      
      # Check if weight data is compatible with country data
      w_pos <- which(as.numeric(dimnames(W[[j]])[[3]]) %in% stats::time(country_data[[j]]))
      x <- matrix(NA, dim(W[[j]])[1], t)
      for (i in 1:t){
        x[, i] <- W[[j]][,, w_pos[i]] %*% X_temp[, i]
      } 
    }
    dimnames(x)[[1]] <- dimnames(W[[j]])[[1]]
    Z <- c(Z, list(x))
    rm(x)
  }
  names(Z) <- countries
  rm(X_temp)
  
  #### Generate data for each country ----
  X_index <- stats::tsp(X)
  data <- c()
  for (i in 1:length(Z)){
    # Split z into domestic and foreign
    x <- stats::ts(t(matrix(Z[[i]][model_specs[[i]]$domestic$variables,],
                            length(model_specs[[i]]$domestic$variables))),
                   start = X_index[1], frequency = X_index[3], class = c("mts", "ts", "matrix"))
    dimnames(x)[[2]] <- model_specs[[i]]$domestic$variables
    
    x_star <- stats::ts(t(matrix(Z[[i]][-(1:length(model_specs[[i]]$domestic$variables)),],
                                 length(model_specs[[i]]$foreign$variables))),
                        start = X_index[1], frequency = X_index[3], class = c("mts", "ts", "matrix"))
    dimnames(x_star)[[2]] <- model_specs[[i]]$foreign$variables
    
    # Check for proper rank specifications
    if (model_specs[[i]]$type == "VEC" & !any(is.na(model_specs[[i]]$cointegration$rank))) {
      k_temp <- length(model_specs[[i]]$domestic$variables)
      if (length(model_specs[[i]]$cointegration$rank) > (k_temp + 1)) {
        model_specs[[i]]$cointegration$rank <- 0:k_temp
      }
      rm(k_temp)
    }
    
    # Collect data
    temp <- list(data = list(domestic = x,
                             foreign = x_star),
                 model = model_specs[[i]])
    
    # Add global data to data element
    if (!is.null(model_specs[[i]]$global)){
      temp$data$global <- X_global
    }
    
    # Add data on deterministic terms ----
    if (!is.null(model_specs[[i]]$deterministic)) {
      
      tt <- nrow(x)
      
      # For VAR models
      if (model_specs[[i]]$type == "VAR") {
        determ <- NULL
        determ_names <- NULL
        if ("const" %in% model_specs[[i]]$deterministic) {
          determ <- cbind(determ, rep(1, tt))
          determ_names <- append(determ_names, "const")
        }
        if ("trend" %in% model_specs[[i]]$deterministic) {
          determ <- cbind(determ, 1:tt)
          determ_names <- append(determ_names, "trend")
        }
        if ("seasonal" %in% model_specs[[i]]$deterministic) {
          freq <- stats::frequency(temp$data$domestic)
          if (freq == 1) {
            warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
          } else {
            # Find first reference date
            pos <- which(stats::cycle(temp$data$domestic) == 1)[1]
            # Define positions which get a seasonal dummy
            pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
            # Produce dummy series
            for (j in 1:(freq - 1)) {
              s_temp <- rep(0, freq)
              s_temp[pos[j]] <- 1
              determ <- cbind(determ, rep(s_temp, length.out = tt))
              determ_names <- c(determ_names, paste("season.", j, sep = ""))
            }
          }
        }
        
        # Make det series a ts object
        determ <- stats::ts(as.matrix(determ), class = c("mts", "ts", "matrix"))
        stats::tsp(determ) <- stats::tsp(x)
        dimnames(determ)[[2]] <- determ_names
      }
      
      # For VEC models
      if (model_specs[[i]]$type == "VEC") {
        
        temp_res <- NULL
        temp_unres <- NULL
        temp_names_res <- NULL
        temp_names_unres <- NULL
        
        use_res <- !is.null(model_specs[[i]]$deterministic$restricted)
        use_unres <- !is.null(model_specs[[i]]$deterministic$unrestricted)
        
        if (use_res) {
          if ("const" %in% model_specs[[i]]$deterministic$restricted) {
            temp_res <- cbind(temp_res, rep(1, tt))
            temp_names_res <- append(temp_names_res, "const")
          }
          if ("trend" %in% model_specs[[i]]$deterministic$restricted) {
            temp_res <- cbind(temp_res, 1:tt)
            temp_names_res <- append(temp_names_res, "trend")
          }
          if ("seasonal" %in% model_specs[[i]]$deterministic$restricted) {
            seas <- NULL
            freq <- stats::frequency(temp$data$domestic)
            if (freq == 1) {
              warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
            } else {
              # Find first reference date
              pos <- which(stats::cycle(temp$data$domestic) == 1)[1]
              # Define positions which get a seasonal dummy
              pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
              for (j in 1:(freq - 1)) {
                s_temp <- rep(0, freq)
                s_temp[pos[j]] <- 1
                temp_res <- cbind(temp_res, rep(s_temp, length.out = tt))
                temp_names_res <- c(temp_names_res, paste("season.", j, sep = "")) 
              }
            }
          }
        }
        
        if (use_unres) {
          if ("const" %in% model_specs[[i]]$deterministic$unrestricted) {
            temp_unres <- cbind(temp_unres, rep(1, tt))
            temp_names_unres <- append(temp_names_unres, "const")
          }
          if ("trend" %in% model_specs[[i]]$deterministic$unrestricted) {
            temp_unres <- cbind(temp_unres, 1:tt)
            temp_names_unres <- append(temp_names_unres, "trend")
          }
          if ("seasonal" %in% model_specs[[i]]$deterministic$unrestricted) {
            seas <- NULL
            freq <- stats::frequency(temp$data$domestic)
            if (freq == 1) {
              warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
            } else {
              # Find first reference date
              pos <- which(stats::cycle(temp$data$domestic) == 1)[1]
              # Define positions which get a seasonal dummy
              pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
              for (j in 1:(freq - 1)) {
                s_temp <- rep(0, freq)
                s_temp[pos[j]] <- 1
                temp_unres <- cbind(temp_unres, rep(s_temp, length.out = tt))
                temp_names_unres <- c(temp_names_unres, paste("season.", j, sep = "")) 
              }
            }
          }
        }
        
        # Make det series ts objects
        determ <- list()
        if (use_res) {
          temp_res <- stats::ts(as.matrix(temp_res), class = c("mts", "ts", "matrix"))
          stats::tsp(temp_res) <- stats::tsp(x)
          dimnames(temp_res)[[2]] <- temp_names_res
          determ$restricted <- temp_res
          rm(temp_res)
        }
        if (use_unres) {
          temp_unres <- stats::ts(as.matrix(temp_unres), class = c("mts", "ts", "matrix"))
          stats::tsp(temp_unres) <- stats::tsp(x)
          dimnames(temp_unres)[[2]] <- temp_names_unres
          determ$unrestricted <- temp_unres
          rm(temp_unres)
        }
      }
      
      temp$data$deterministic <- determ
      
    }
    
    temp$data$weights = W[[i]]
    
    data <- c(data, list(temp))
    rm(list = c("x", "x_star"))
  }
  names(data) <- countries
  
  #### Generate country model for each lag and rank specification ####
  data_temp <- NULL
  names_temp <- NULL
  for (i in countries){
    
    # Helping object, since loops cannot deal with empty objects
    if (!global) {
      data[[i]]$model$global$lags <- "none"
    }
    if (data[[i]]$model$type == "VAR") {
      data[[i]]$model$rank <- "none"
    }
    
    # Create country objects
    for (j in data[[i]]$model$domestic$lags) {
      for (k in data[[i]]$model$foreign$lags) {
        for (l in data[[i]]$model$global$lags) {
          for (m in data[[i]]$model$rank) {
            
            temp <- data[[i]]
            temp$model$domestic$lags <- j
            temp$model$foreign$lags <- k
            if (global) {
              temp$model$global$lags <- l
            } else {
              temp$model$global <- NULL
            }
            
            if (temp$model[["type"]] == "VEC") {
              temp$model[["rank"]] <- m
            } else {
              temp$model[["rank"]] <- NULL
            }
            
            temp$model$varselect <- "none"
            
            # Add matrices for estimation
            if (temp$model$type == "VAR") {
              temp_i <- .gen_varx(temp)
              temp$data[["Y"]] <- temp_i[["Y"]]
              temp$data[["Z"]] <- temp_i[["Z"]]
              temp$data[["SUR"]] <- temp_i[["SUR"]]
            }
            if (temp$model$type == "VEC") {
              temp_i <- .gen_vecx(temp)
              temp$data[["Y"]] <- temp_i[["Y"]]
              temp$data[["W"]] <- temp_i[["W"]]
              temp$data[["X"]] <- temp_i[["X"]]
              temp$data[["SUR"]] <- temp_i[["SUR"]]
            }
            
            data_temp <- c(data_temp, list(temp))
            names_temp <- c(names_temp, i)
          } 
        }
      }
    }
  }
  
  names(data_temp) <- names_temp
  
  type <- unique(unlist(lapply(data_temp, function(x){x$model$type})))
  if (length(type) > 1) {
    warning("Multiple types of submodels identified.")
  }
  if (type[1] == "VAR") {
    class(data_temp) <- append("gvarsubmodels", class(data_temp))  
  }
  if (type[1] == "VEC") {
    class(data_temp) <- append("gvecsubmodels", class(data_temp))  
  }
  
  return(data_temp)
}
