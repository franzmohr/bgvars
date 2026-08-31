#' Create Submodel Specifications
#' 
#' Produces a list of model specifications for each entity in a GVAR model.
#' 
#' @param object an object of class 'gvarmodel' or 'gvecmodel'.
#' @param submodel name of the sub-model, for which input data should be generated
#' based on the input in argument \code{object}.
#' @param endogen character vector of variables that should enter each sub-model
#' as endogenous variables, if they are available for the respective sub-model.
#' @param p_endogen an integer vector of the lag order (default is \code{p_endogen = 1})
#' of a sub-model's endogenous variables.
#' @param exogen character vector of variables that should enter each sub-model
#' as weakly exogenous variables.
#' @param p_exogen an integer vector of the lag order (default is \code{p_exogen = 1})
#' of a sub-model's weakly exogenous variables.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the lag order of a sub-model's global variables.
#' If \code{NULL} (default), models do not include global variables.
#' @param deterministic a character specifying which deterministic terms should
#' be included. Available values are \code{"none"}, \code{"const"} (default) for an intercept,
#' \code{"trend"} for a linear trend, and \code{"both"} for an intercept with a linear trend.
#' @param seasonal logical. If \code{TRUE}, seasonal dummy variables are
#' generated as additional deterministic terms. The amount of dummies depends on the frequency of the
#' time-series object provided in \code{object}. Defaults to \code{FALSE}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model. Defaults to \code{FALSE}.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details
#' 
#' Argument \code{error} specifies the structure of the covariance matrix of
#' the error term and how it is estimated. Possible specifications are:
#' \itemize{
#'  \item{\code{"wishart"}: The covariance is estimated using a Wishart prior.}
#'  \item{\code{"gamma"}: Only the diagonal elements of the covariance matrix
#'  are estimated using a gamma prior.
#'  Off-diagonal elements are not estimated and set to zero.}
#'  \item{\code{"gamma+covar"}: The diagonal elements of the covariance matrix
#'  are estimated using a gamma prior.
#' Covariances are estimated based on a triangular decomposition.}
#'  \item{\code{"sv"}: Only the diagonal elements of the covariance matrix are
#'  estimated using a stochastic volatility
#' algorithm. Off-diagonal elements are not estimated and set to zero.}
#'  \item{\code{"sv+covar"}: Only the diagonal elements of the covariance matrix
#'  are estimated using a stochastic volatility
#' algorithm. Covariances are estimated based on a triangular decomposition.}
#' }
#' 
#' Available specifications for argument \code{varsel} are:
#' \itemize{
#'  \item{\code{"none"}: No variable selection algorithm is used.}
#'  \item{\code{"bvs"}: Bayesian variable selection as proposed in Korobilis (2013).}
#'  \item{\code{"ssvs"}: Stochastic search variable selection as proposed in George et al. (2008).}
#' }
#' 
#' @return An object of class 'modellist'.
#' 
#' 
#' @examples
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]]
#' global_data <- gvar2023[["global_data"]]
#' 
#' # Create empty model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#' 
#' # Create sub-models
#' model <- create_varxsubmodel(object,
#'                              submodel = "AT",
#'                              endogen = c("y","p", "rs"), p_endogen = 1,
#'                              exogen = c("y", "p"), p_exogen = 1,
#'                              global = "poil", s = 0,
#'                              deterministic = "const", seasonal = FALSE,
#'                              structural = FALSE, tvp = FALSE,
#'                              error = "wishart", varsel = "none",
#'                              iterations = 10000, burnin = 2000)
#' 
#' 
#' @export
create_varxsubmodel <- function(object,
                                submodel,
                                endogen = NULL,
                                p_endogen = 1,
                                exogen = NULL,
                                p_exogen = 1,
                                global = NULL,
                                s = NULL,
                                deterministic = "const",
                                seasonal = FALSE,
                                structural = FALSE,
                                tvp = FALSE,
                                error = "wishart",
                                varsel = "none",
                                iterations = 10000,
                                burnin = 2000){
  
  if (length(submodel) > 1) {
    stop("Argument 'submodel' may only contain one element.")
  }
  
  tt <- nrow(object[["global"]][["endogen"]])
  tsp_global <- stats::tsp(object[["global"]][["endogen"]])
  index <- object[["global"]][["index"]]
  
  # Get n_z
  vars_endogen_old <- index[which(index[, "submodel"] == submodel) , "variable"]
  n_endogen_old <- length(vars_endogen_old)
  vars_exogen_old <- unique(index[index[, "submodel"] != submodel, "variable"])
  n_exogen_old <- length(vars_exogen_old)
  n_z_old <- n_endogen_old + n_exogen_old
  
  # Get positions for submodel
  if (is.null(endogen)) {
    pos_endogen <- which(index[, "submodel"] == submodel)
    vars_endogen <- index[pos_endogen, "variable"]
  } else {
    pos_endogen <- which(index[, "submodel"] == submodel & index[, "variable"] %in% endogen)
    if (length(pos_endogen) == 0) {
      stop(paste0("For submodel ", submodel, " no variable from argument 'endogen' is available."))
    }
    vars_endogen <- index[pos_endogen, "variable"]
  }
  pos_endogen <- which(vars_endogen_old %in% vars_endogen)
  n_endogen <- length(vars_endogen)
  
  
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])
  if (is.null(exogen)) {
    pos_exogen <- n_endogen_old + 1:length(vars_exogen)
    n_exogen <- length(vars_exogen)
  } else {
    pos_exogen <- which(vars_exogen %in% exogen)
    vars_exogen <- vars_exogen[pos_exogen]
    n_exogen <- length(vars_exogen)
    pos_exogen <- n_endogen_old + pos_exogen
  }
  n_z <- n_endogen + n_exogen
  
  pos_new <- c(pos_endogen, pos_exogen)
  
  # Create z
  z_i <- matrix(NA, n_z, tt)
  for (i in 1:tt) {
    z_i[, i] <- object[["weights"]][[submodel]][(i - 1) * n_z_old + pos_new,] %*% object[["global"]][["endogen"]][i, ]
  }
  z_i <- stats::ts(t(z_i))
  stats::tsp(z_i) <- tsp_global
  
  # Endogenous variables
  endogen <- z_i[, 1:n_endogen]
  if (is.null(dimnames(endogen))) {
    # If 'endogen' is a simple ts object, transform it into a matrix object
    # to keep variable name information
    endogen <- stats::ts(as.matrix(endogen), class = c("mts", "ts", "matrix"))
    stats::tsp(endogen) <- tsp_global
  }
  dimnames(endogen)[[2]] <- vars_endogen
  
  
  # Weakly exogenous variables
  exogen <- z_i[, n_endogen + 1:n_exogen]
  if (is.null(dimnames(exogen))) {
    # If 'exogen' is a simple ts object, transform it into a matrix object
    # to keep variable name information
    exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
    stats::tsp(exogen) <- tsp_global
  }
  vars_exogen <- paste0(vars_exogen, ".s")
  dimnames(exogen)[[2]] <- vars_exogen
  
  # Global variables
  use_global <- !is.null(object[["global"]][["exogen"]]) & !is.null(global)
  n_global <- 0L
  if (use_global) {
    if (is.null(s)) {
      stop("If global variables should be used, argument 's' must be specified.") 
    }
    pos_global <- which(dimnames(object[["global"]][["exogen"]])[[2]] %in% global)
    vars_global <- dimnames(object[["global"]][["exogen"]])[[2]][pos_global]
    global <- object[["global"]][["exogen"]][, pos_global]
    n_global <- length(vars_global)
    
    if (is.null(dimnames(global))) {
      # If 'global' is a simple ts object, transform it into a matrix object
      # to keep variable name information
      tsp_temp <- stats::tsp(global)
      global <- stats::ts(as.matrix(global), class = c("mts", "ts", "matrix"))
      stats::tsp(global) <- tsp_temp
    }
    dimnames(global)[[2]] <- vars_global
  }
  
  # From here on its nearly identicl to create_bvarmodel
  if (n_endogen == 1 & structural) {
    # Overwrite structural parameter if there is only one endogenous variable
    structural <- FALSE
    if (error == "gamma+covar") {
      error <- "gamma"
    }
    if (error == "sv+covar") {
      error <- "sv"
    }
  }
  
  
  if (structural & error %in% c("wishart", "gamma+covar", "sv+covar")) {
    stop(paste0("Structural models cannot be estimated with argument 'error' specified as '", error,"'."))
  }
  
  if (!varsel %in% c("none", "bvs", "ssvs")) {
    stop("Specification of argument 'varsel' is not supported.")
  }
  
  algo <- NULL
  if (tvp) {
    algo <- paste0(algo, "Tvp")
  } else {
    algo <- paste0(algo, "Normal")
  }
  if (error == "wishart") {
    algo <- paste0(algo, "Wishart")
  }
  if (error %in% c("gamma", "gamma+covar")) {
    algo <- paste0(algo, "Gamma")
  }
  if (error %in% c("sv", "sv+covar")) {
    algo <- paste0(algo, "Stochvol")
  }
  algo <- paste0("Var", algo)
  
  model <- NULL
  model[["type"]] <- ifelse(n_endogen == 1, "ARX", "VARX")
  model[["algorithm"]] <- algo
  model[["k"]] <- n_endogen
  model[["p"]] <- 0L
  model[["m"]] <- 0L
  model[["s"]] <- 0L
  model[["n"]] <- 0L
  model[["k_endogen"]] <- n_endogen
  model[["p_endogen"]] <- 0L
  model[["k_exogen"]] <- n_exogen
  model[["p_exogen"]] <- 0L
  model[["m_global"]] <- n_global
  model[["s_global"]] <- 0L
  model[["varsel"]] <- varsel
  model[["endogen"]] <- vars_endogen
  model[["exogen"]] <- vars_exogen
  if (use_global) {
    model[["global"]] <- vars_global 
  }
  
  
  p_endogen_max <- max(p_endogen)
  temp <- endogen
  temp_name <- vars_endogen
  if (p_endogen_max >= 1) {
    # Obtain lags of endogenous variables
    for (i in 1:p_endogen_max) {
      temp <- cbind(temp, stats::lag(endogen, -i))
      if (nchar(p_endogen_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_endogen_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0(vars_endogen, ".", i_temp))
    }
  }
  
  
  p_exogen_max <- max(p_exogen)
  temp <- cbind(temp, exogen)
  if (nchar(p_exogen_max) > 2) {
    i_temp <- rep(0, nchar(p_exogen_max))
  } else {
    i_temp <- rep(0, 2)
  }
  i_temp <- paste0(i_temp, collapse = "")
  temp_name <- c(temp_name, paste0(vars_exogen, ".l", i_temp))
  if (p_exogen_max > 0) {
    for (i in 1:p_exogen_max) {
      temp <- cbind(temp, stats::lag(exogen, -i))
      if (nchar(p_exogen_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_exogen_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0(vars_exogen, ".l", i_temp))
    } 
  }
  
  if (use_global) {
    s_max <- max(s)
    temp <- cbind(temp, global)
    if (nchar(s_max) > 2) {
      i_temp <- rep(0, nchar(s_max))
    } else {
      i_temp <- rep(0, 2)
    }
    i_temp <- paste0(i_temp, collapse = "")
    temp_name <- c(temp_name, paste0(vars_global, ".l", i_temp))
    if (s_max > 0) {
      for (i in 1:s_max) {
        temp <- cbind(temp, stats::lag(global, -i))
        if (nchar(s_max) > 2) {
          i_temp <- paste0(c(rep(0, nchar(s_max) - nchar(i)), i), collapse = "")
        } else {
          i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
        }
        temp_name <- c(temp_name, paste0(vars_global, ".l", i_temp))
      } 
    }
  } else {
    s_max <- 0
  }
  
  # Determinsitic terms
  det_data <- NULL
  det_name <- NULL
  det_pos <- ncol(temp)
  
  # Add intercept term
  if (deterministic %in% c("const", "both")) {
    temp <- cbind(temp, 1)
    temp_name <- c(temp_name, "const")
    det_name <- c(det_name, "const")
  }
  
  # Add linear trend
  if (deterministic %in% c("trend", "both")) {
    temp <- cbind(temp, 1:nrow(temp) - max(p_endogen_max, p_exogen_max, s_max))
    temp_name <- c(temp_name, "trend")
    det_name <- c(det_name, "trend")
  }
  
  # Add seasonal dummies
  if (seasonal) {
    freq <- stats::frequency(endogen)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(stats::cycle(temp) == 1)[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      for (i in 1:(freq - 1)) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = nrow(temp)))
        temp_name <- c(temp_name, paste("season.", i, sep = ""))
        det_name <- c(det_name, paste("season.", i, sep = ""))
      }
    }
  }
  
  # Update model specs for deterministic terms
  use_det <- FALSE
  if (length(det_name) > 0) {
    model[["n"]] <- length(det_name)
    use_det <- TRUE
    det_data <- temp[, det_pos + 1:model[["n"]]]
    
    # If 'det_data' is a simple ts object, transform it into a matrix object
    # to keep variable name information
    tsp_det_data <- stats::tsp(det_data)
    det_data <- stats::ts(as.matrix(det_data), class = c("mts", "ts", "matrix"))
    stats::tsp(det_data) <- tsp_det_data
    dimnames(det_data) <- list(NULL, det_name)
  }
  
  temp <- stats::na.omit(temp)
  
  # Set if the model is structural
  if ("logical" %in% class(structural)) {
    model[["structural"]] <- structural
    if (structural) {
      model[["type"]] <- "SVARX" 
    }
  } else {
    stop("Argument 'structural' must be of class 'logical'.")
  }
  
  ## errors ----
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      stop("Invalid specification of argument 'error'.")
    }
    model[["error"]] <- error
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  ## tvp ----
  if ("logical" %in% class(tvp)) {
    model[["tvp"]] <- tvp
  } else {
    stop("Argument 'tvp' must be of class 'logical'.")
  }
  
  # Iterations and burnin ----
  model[["iterations"]] <- as.integer(iterations)
  model[["burnin"]] <- as.integer(burnin)
  
  # Endogenous variables y
  y <- stats::ts(as.matrix(temp[, 1:n_endogen]), class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- stats::tsp(temp)
  dimnames(y)[[2]] <- temp_name[1:n_endogen]
  
  # Structural data
  y_A0 <- NULL
  if (structural & n_endogen > 1) {
    y_A0 <- kronecker(-y, diag(1, n_endogen))
    pos <- NULL
    for (j in 1:n_endogen) {
      pos <- c(pos, (j - 1) * n_endogen + 1:j)
    }
    y_A0 <- y_A0[, -pos]
  }
  
  # Use fake loop iteration range to handle modesl without global variables
  if (!use_global) {
    s <- 99
  }
  
  result <- NULL
  for (i in p_endogen) { # for each lag p_endogen
    for (j in p_exogen) { # for each lag p_exogen
      for (k in s) {
        pos <- NULL
        model_i <- model
        if (i >= 1) {
          pos <- c(pos, n_endogen + 1:(n_endogen * i))
          model_i[["p"]] <- as.integer(i)
          model_i[["p_endogen"]] <- as.integer(i)
        }  
        
        pos <- c(pos, n_endogen + n_endogen * p_endogen_max + 1:(n_exogen * (j + 1)))
        model_i[["p_exogen"]] <- as.integer(j)
        
        if (use_global) {
          pos <- c(pos, n_endogen + n_endogen * p_endogen_max + n_exogen * (p_exogen_max + 1) + 1:(n_global * (k + 1)))
          model_i[["m_global"]] <- as.integer(n_global)
          model_i[["s_global"]] <- as.integer(k)
        }
        
        model_i[["m"]] <- model_i[["k_exogen"]] * (model_i[["p_exogen"]] + 1) + model_i[["m_global"]] * (model_i[["s_global"]] + 1)
        model_i[["s"]] <- 0L
        
        if (use_det) {
          pos <- c(pos, n_endogen + n_endogen * p_endogen_max + n_exogen * (p_exogen_max + 1) + n_global * (s_max + 1) + 1:length(det_name))
        }
        
        x <- NULL
        z <- NULL
        if (length(pos) > 0) {
          # Create data input matrix of the respective model
          x <- stats::ts(as.matrix(temp[, pos]), class = c("mts", "ts", "matrix")) 
          stats::tsp(x) <- stats::tsp(temp)
          dimnames(x)[[2]] <- temp_name[pos]
          z <- kronecker(x, diag(1, n_endogen))
        }
        
        # If specified add structural data to SUR form
        if (!is.null(y_A0)) {
          z <- cbind(z, y_A0)
        }
        dimnames(z) <- NULL
        
        # Create individual model
        result_i <- list("model" = model_i,
                         "data" = list("original" = list("endogen" = endogen,
                                                         "exogen" = exogen,
                                                         "deterministic" = det_data),
                                       "train" = list("y" = y,
                                                      "x" = x,
                                                      "z" = z)))
        
        # Update class of individual model
        class(result_i) <- c("varxsubmodel", "bvarmodel", "list")
        
        result <- c(result, list(result_i)) 
      }
    }
  }
  
  class(result) <- c("modellist", "list")
  
  return(result)
}