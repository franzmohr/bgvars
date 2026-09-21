#' Create Sub-Models
#' 
#' Produces a list of VECX models for each sub-model in a GVEC model.
#' 
#' @param object an object of class 'gvarmodel' or 'gvecmodel'.
#' @param submodel name of the sub-model, for which input data should be generated
#' based on the input in argument \code{object}.
#' @param endogen character vector of variables that should enter each sub-model
#' as endogenous variables, if they are available.
#' @param p_endogen an integer vector of the lag order of a sub-model's
#' endogenous variables in the corresponding VAR model in levels (default is
#' \code{p_endogen = 1}). It must be at least 1, since the sub-model holds
#' \code{p_endogen - 1} lagged differences of the endogenous variables.
#' @param exogen character vector of variables that should enter each sub-model
#' as weakly exogenous variables.
#' @param p_exogen an integer vector of the lag order of a sub-model's weakly
#' exogenous variables in the corresponding VARX sub-model in levels, counted
#' from the contemporaneous term onwards as in \code{\link{create_varxsubmodel}}
#' (default is \code{p_exogen = 2}). It must be at least 1. The lagged levels
#' \eqn{x_{t-1}} enter the error correction term whatever the value, and the
#' sub-model holds the \code{p_exogen - 1} differences from \eqn{\Delta x_t}
#' onwards, so 1 uses the lagged levels alone, 2 adds \eqn{\Delta x_t} and 3
#' adds \eqn{\Delta x_{t-1}}. A VARX\eqn{^*}(p, q) model, whose weakly exogenous
#' variables reach back to lag q, is therefore specified by
#' \code{p_endogen = p} and \code{p_exogen = q + 1}, and
#' \code{\link[bvartools]{vec_to_var}} returns a VARX sub-model with the same
#' \code{p_exogen}. The exception is \code{p_exogen = 1}, whose lagged levels
#' take two blocks in levels, the contemporaneous one with a coefficient of zero.
#' @param global character vector of variables that should enter each sub-model
#' as global variables.
#' @param s an integer vector of the lag order of a sub-model's global
#' variables in the corresponding VARX sub-model in levels, counted in the same
#' way as \code{p_exogen}. It must be at least 1. If \code{NULL} (default),
#' models do not include global variables.
#' @param r an integer vector of the cointegration rank. See 'Details'.
#' @param const a character specifying whether a constant term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param trend a character specifying whether a trend term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param seasonal a character specifying whether seasonal dummies should be included in the error
#' correction term (\code{"restricted"}) or in the non-cointegreation term as \code{"unrestricted"}
#' variables. If \code{NULL} (default) no seasonal terms will be added. The amount of dummy variables
#' will be automatically detected and depends on the frequency of the time-series object provided
#' in \code{data}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VEC model. Defaults to \code{FALSE}.
#' @param tvp logical indicating whether the model parameters are time varying.
#' Defaults to \code{FALSE}.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param delta_beta,delta_sigma numeric discount factors in \eqn{(0, 1]} of the
#' discounted model, the first governing the coefficients and the second the
#' error covariance. Both default to 1, at which the quantity they govern does
#' not move. Ignored unless \code{algorithm = "discount"}, and a vector in
#' either produces one sub-model per value. See 'Details'.
#' @param algorithm algorithm that should be used for posterior simulation. If \code{NULL}
#' (default), standard algorithms will be used. See 'Details' for available
#' non-standard options.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 10000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' @param thin an integer thinning interval of the sampler (defaults to 1). After
#' the burn-in the sampler keeps the last of every \code{thin} draws, so it runs
#' \code{burnin + iterations * thin} draws and still keeps \code{iterations}.
#' Unlike \code{\link[coda]{thin}}, which thins draws already made, the
#' draws that are not kept are never held in memory, which is what lets long
#' chains of time varying sub-models fit into it.
#'
#' @details
#' 
#' If an integer vector is provided as argument \code{p_endogen}, \code{p_exogen},
#' \code{s} or \code{r}, the function will produce a distinct model for all
#' possible combinations of those specifications.
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
#' Available specifications for argument \code{algorithm} are:
#' \itemize{
#'    \item{\code{"KLGS2010"}: Algorithm proposed in Koop, León-González & Strachan (2010).}
#'    \item{\code{"discount"}: The discounted time varying parameter model of
#'    West & Harrison (1997, ch. 16) with the discounted Wishart of Uhlig
#'    (1997), estimated by \code{VecTvpDiscount}.}
#' }
#'
#' The discounted model is not a sampler. Its posterior is closed form -- one
#' pass over the sample, no chain -- so \code{burnin} must be 0 and \code{thin}
#' 1, \code{error} must be \code{"wishart"}, and neither variable selection nor
#' a structural model is available. It conditions on a fixed cointegration
#' matrix, which \code{\link[bvartools]{add_initial_values}} takes from
#' Johansen's estimate, so what drifts is the adjustment to the long-run
#' relations and not the relations themselves. What it buys is speed and an
#' exact marginal likelihood: the sum of \code{/posterior/loglik} is the log
#' marginal likelihood of the sub-model given its rank, its cointegration
#' matrix, its lag orders and the two discounts, so a grid over any of them is
#' compared without a chain being run for any of it.
#'
#' A sub-model of that algorithm is estimated by BayesTS on its file, as every
#' sub-model of an export is, and by nothing else: there is no R implementation
#' of the filter. \code{\link[bvartools]{is_discount_model}} says which
#' algorithms these are.
#' 
#' @return An object of class 'modellist', which contains at least one element
#' of class 'vecxsubmodel'.
#' 
#' 
#' @examples
#' 
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("EA", "CA", "US"))
#' 
#' # Create empty model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#' 
#' # Create sub-models
#' model <- create_vecxsubmodel(object,
#'                              submodel = "EA",
#'                              r = 1,
#'                              global = "poil",
#'                              s = 2,
#'                              iterations = 10,
#'                              burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' 
#' @export
create_vecxsubmodel <- function(object,
                                submodel,
                                endogen = NULL,
                                p_endogen = 1,
                                exogen = NULL,
                                p_exogen = 2,
                                global = NULL,
                                s = NULL,
                                r = NULL,
                                const = NULL,
                                trend = NULL,
                                seasonal = NULL,
                                structural = FALSE,
                                tvp = FALSE,
                                error = "wishart",
                                varsel = "none",
                                algorithm = NULL,
                                delta_beta = 1,
                                delta_sigma = 1,
                                iterations = 10000,
                                burnin = 2000,
                                thin = 1){

  if (length(submodel) > 1) {
    stop("Argument 'submodel' may only contain one element.")
  }
  
  if (any(p_endogen < 1)) {
    stop("Argument 'p_endogen' must be at least 1 for VECX models.")
  }
  
  if (any(p_exogen < 1)) {
    stop("Argument 'p_exogen' must be at least 1 for VECX models.")
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
  
  # Get positions of sub-model's variables in its weight matrix
  vars_endogen <- index[index[, "submodel"] == submodel, "variable"]
  if (is.null(endogen)) {
    pos_endogen <- 1:length(vars_endogen)
  } else {
    # Endogenous variables are used if they are available, so the ones that are
    # not are dropped rather than objected to. match() answers with NA for
    # those, and testing the length of its result would never find them: it is
    # the length of 'endogen'.
    pos_endogen <- match(endogen, vars_endogen)
    pos_endogen <- pos_endogen[!is.na(pos_endogen)]
    if (length(pos_endogen) == 0) {
      stop(paste0("For sub-model ", submodel, " no variable from argument 'endogen' is available."))
    }
    vars_endogen <- vars_endogen[pos_endogen]
  }
  n_endogen <- length(vars_endogen)
  
  
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])
  if (is.null(exogen)) {
    pos_exogen <- n_endogen_old + 1:length(vars_exogen)
  } else {
    pos_exogen <- match(exogen, vars_exogen)
    if (any(is.na(pos_exogen))) {
      stop(paste0("For sub-model ", submodel, " at least one specified exogenous variable is not available."))
    }
    vars_exogen <- vars_exogen[pos_exogen]
    pos_exogen <- n_endogen_old + pos_exogen
  }
  n_exogen <- length(vars_exogen)
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
    if (any(s < 1)) {
      stop("Argument 's' must be at least 1 for VECX models.")
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
  algo <- paste0("Vec", algo)

  # The one algorithm outside the naming grammar that a sub-model can use. The
  # specification is refused by bvartools, which raises what BayesTS would raise
  # against the written file -- and raises it once for a grid of sub-models
  # rather than once per file.
  if (!is.null(algorithm)) {
    if (!identical(algorithm, "discount")) {
      stop("Specified algorithm not recognized.")
    }
    algo <- "VecTvpDiscount"
    bvartools::check_discount_specification(
      k = n_endogen, error = error, varsel = varsel, structural = structural,
      burnin = burnin, thin = thin, delta_beta = delta_beta,
      delta_sigma = delta_sigma)
    # The coefficients of this model drift, which is what `tvp` says. How far
    # they drift is what `delta_beta` says, and it is a model in its own right
    # at one, where they do not move at all.
    tvp <- TRUE
  }
  use_discount <- identical(algo, "VecTvpDiscount")
  
  model <- NULL
  model[["type"]] <- ifelse(n_endogen == 1, "ECX", "VECX")
  model[["algorithm"]] <- algo
  model[["k"]] <- n_endogen
  model[["p"]] <- 0L
  model[["m"]] <- 0L
  model[["s"]] <- 0L
  model[["n"]] <- 0L
  model[["n_restricted"]] <- 0L
  model[["rank"]] <- 0L
  model[["k_beta"]] <- 0L
  
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
  
  # ****************************************************************************
  # Sample creation ----
  
  # Differenced endogenous variables
  diff_endogen <- diff(endogen)
  temp_name <- paste("d.", vars_endogen, sep = "")
  temp <- diff_endogen
  
  # Endogenous ECT variables
  temp <- cbind(temp, stats::lag(endogen, -1))
  temp_name <- c(temp_name, paste("l.", vars_endogen, sep = ""))
  n_ect <- length(vars_endogen)
  
  # Weakly exogenous ECT variables
  temp <- cbind(temp, stats::lag(exogen, -1))
  temp_name <- c(temp_name, paste("l.", vars_exogen, sep = ""))
  n_ect <- n_ect + length(vars_exogen)
  
  # Global ECT variables
  if (use_global) {
    temp <- cbind(temp, stats::lag(global, -1))
    temp_name <- c(temp_name, paste("l.", vars_global, sep = ""))
    n_ect <- n_ect + length(vars_global)
  }
  
  # Lags of differenced endogenous variables
  p_endogen_max <- max(p_endogen)
  if (p_endogen_max > 1) {
    # Obtain lags of endogenous variables
    for (i in 1:(p_endogen_max - 1)) {
      temp <- cbind(temp, stats::lag(diff_endogen, -i))
      if (nchar(p_endogen_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_endogen_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0("d.",vars_endogen, ".", i_temp))
    }
  }
  
  # Lags of differenced weakly exogenous variables
  #
  # 'p_exogen' is the lag order of the VARX sub-model in levels, counted from
  # the contemporaneous term onwards, so the VECX sub-model holds one block of
  # differences less: p_exogen - 1 of them, from the contemporaneous one on.
  d_exogen_max <- max(p_exogen) - 1
  diff_exogen <- diff(exogen)
  # A lag order of one asks for no differences of the weakly exogenous
  # variables at all, so the contemporaneous block must not be added either.
  # The offsets used below count this block as n_exogen * d_exogen_max columns,
  # and adding it unconditionally made them address the wrong columns.
  if (d_exogen_max > 0) {
    temp <- cbind(temp, diff_exogen)
    if (nchar(d_exogen_max) > 2) {
      i_temp <- rep(0, nchar(d_exogen_max))
    } else {
      i_temp <- rep(0, 2)
    }
    i_temp <- paste0(i_temp, collapse = "")
    temp_name <- c(temp_name, paste0("d.", vars_exogen, ".", i_temp))
  }
  if (d_exogen_max > 1) {
    for (i in 1:(d_exogen_max - 1)) {
      temp <- cbind(temp, stats::lag(diff_exogen, -i))
      if (nchar(d_exogen_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(d_exogen_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0("d.", vars_exogen, ".", i_temp))
    } 
  }
  
  if (use_global) {
    d_global_max <- max(s) - 1
    diff_global <- diff(global)
    # As for the weakly exogenous variables, s = 1 means that the global
    # variables enter the error correction term in levels but contribute no
    # differences to the non-cointegration term.
    if (d_global_max > 0) {
      temp <- cbind(temp, diff_global)
      if (nchar(d_global_max) > 2) {
        i_temp <- rep(0, nchar(d_global_max))
      } else {
        i_temp <- rep(0, 2)
      }
      i_temp <- paste0(i_temp, collapse = "")
      temp_name <- c(temp_name, paste0("d.", vars_global, ".", i_temp))
    }
    if (d_global_max > 1) {
      for (i in 1:(d_global_max - 1)) {
        temp <- cbind(temp, stats::lag(diff_global, -i))
        if (nchar(d_global_max) > 2) {
          i_temp <- paste0(c(rep(0, nchar(d_global_max) - nchar(i)), i), collapse = "")
        } else {
          i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
        }
        temp_name <- c(temp_name, paste0("d.", vars_global, ".", i_temp))
      } 
    }
  } else {
    d_global_max <- 0
  }
  
  temp <- stats::na.omit(temp)
  ts_info <- stats::tsp(temp)
  
  # Final endogenous variables
  y <- stats::ts(as.matrix(temp[, 1:n_endogen]), class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- ts_info
  dimnames(y)[[2]] <- temp_name[1:n_endogen]
  
  tt <- nrow(temp)
  
  ect <- matrix(temp[, n_endogen + 1:n_ect], tt)
  ect_names <- temp_name[n_endogen + 1:n_ect]
  
  x <- matrix(temp[, -(1:(n_endogen + n_ect))], tt)
  x_names <- temp_name[-(1:(n_endogen + n_ect))]
  
  det_name_r <- NULL
  det_name_ur <- NULL
  n_det_ur <- 0

  # The deterministic terms are taken from the global model rather than built
  # here, so that every sub-model uses the same series whatever sample it ends
  # up with. See .global_deterministic(). The window is the one 'temp' was
  # trimmed to, so the columns line up with 'ect' and 'x' row by row.
  tsp_temp <- stats::tsp(temp)
  det_available <- dimnames(object[["global"]][["deterministic"]])[[2]]
  det_global <- .submodel_deterministic(object, det_available,
                                        start = tsp_temp[1], end = tsp_temp[2])
  # They join 'ect' and 'x' as plain columns. Keeping the time series class
  # would send cbind() to its ts method, which cannot extend a matrix that has
  # no columns yet -- the case of a model whose only regressors are
  # deterministic. The tsp attribute of the results is set further below.
  det_global <- unclass(det_global)
  attr(det_global, "tsp") <- NULL
  
  if (!is.null(const)) {
    if (const == "restricted") {
      ect <- cbind(ect, det_global[, "const"])
      ect_names <- c(ect_names, "const") 
      det_name_r <- c(det_name_r, "const") 
      n_ect <- n_ect + 1
    }
    
    if (const == "unrestricted") {
      x <- cbind(x, det_global[, "const"])
      x_names <- c(x_names, "const")
      det_name_ur <- c(det_name_ur, "const") 
      n_det_ur <- n_det_ur + 1
    }
  }
  
  if (!is.null(trend)) {
    if (trend == "restricted") {
      ect <- cbind(ect, det_global[, "trend"])
      ect_names <- c(ect_names, "trend")
      det_name_r <- c(det_name_r, "trend") 
      n_ect <- n_ect + 1
    }
    
    if (trend == "unrestricted") {
      x <- cbind(x, det_global[, "trend"])
      x_names <- c(x_names, "trend")
      det_name_ur <- c(det_name_ur, "trend")
      n_det_ur <- n_det_ur + 1
    }
  }
  
  if(!is.null(seasonal)) {
    freq <- tsp_global[3]
    if (freq == 1) {
      # Nothing to add, as for a VARX sub-model. The dummies used to be added
      # after the warning all the same, from a variable that was never created.
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      s_name <- paste0("season.", 1:(freq - 1))
      seas <- det_global[, s_name, drop = FALSE]

      if (seasonal == "restricted") {
        ect <- cbind(ect, seas)
        ect_names <- c(ect_names, s_name)
        det_name_r <- c(det_name_r, s_name)
        n_ect <- n_ect + freq - 1
      }

      if (seasonal == "unrestricted") {
        x <- cbind(x, seas)
        x_names <- c(x_names, s_name)
        det_name_ur <- c(det_name_ur, s_name)
        n_det_ur <- n_det_ur + length(s_name)
      }
    }
  }
  
  det_data <- NULL
  
  use_det_ur <- FALSE
  if (length(det_name_ur) > 0) {
    use_det_ur <- TRUE
    model[["n"]] <- length(det_name_ur)
    model[["deterministic"]] <- det_name_ur
    det_data <- x[, which(x_names %in% det_name_ur)]
  }
  
  use_det_r <- FALSE
  if (length(det_name_r) > 0) {
    use_det_r <- TRUE
    model[["n_restricted"]] <- length(det_name_r)
    model[["deterministic_restricted"]] <- det_name_r
    det_data <- cbind(det_data, ect[, which(ect_names %in% det_name_r)])
  }
  
  if (is.null(r)) {
    if (n_ect > n_endogen) {
      r <- 0:n_endogen
    } else {
      r <- 0:(n_endogen - 1)
    }
    message("Argument rank 'r' not specified. Generating models for r = ", paste0(r, collapse = ", "), ".")
  } else {
    if (any(r > n_endogen)) {
      stop("Argument rank 'r' must be smaller than or equal to the number of endogenous variables.")
    }
  }
  
  # Set if the model is structural
  if ("logical" %in% class(structural)) {
    model[["structural"]] <- structural
    if (structural) {
      model[["type"]] <- "SVECX" 
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
  # Carried only when it thins, as bvartools does for its own models, whose
  # samplers and labels of the kept draws read it from here.
  model[["thin"]] <- .check_sampler_thin(thin)
  
  
  ect <- stats::ts(as.matrix(ect), class = c("mts", "ts", "matrix"))
  stats::tsp(ect) <- ts_info
  dimnames(ect)[[2]] <- ect_names
  
  if (length(x_names) > 0) {
    x <- stats::ts(as.matrix(x), class = c("mts", "ts", "matrix"))
    stats::tsp(x) <- ts_info
    dimnames(x)[[2]] <- x_names
  } else {
    x <- NULL
  }
  
  if (!is.null(det_data)) {
    det_data <- stats::ts(as.matrix(det_data), class = c("mts", "ts", "matrix"))
    stats::tsp(det_data) <- ts_info
    dimnames(det_data)[[2]] <- c(det_name_ur, det_name_r)
  }
  
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
  
  # Use fake loop iteration range to handle models without global variables
  if (!use_global) {
    s <- 99
  }
  
  # A grid over the two discounts, for the one algorithm that has them. Every
  # other algorithm gets the single specification it always had, so the two
  # inner loops run once and change nothing.
  if (use_discount) {
    grid_beta <- as.numeric(delta_beta)
    grid_sigma <- as.numeric(delta_sigma)
  } else {
    grid_beta <- NA_real_
    grid_sigma <- NA_real_
  }

  result <- NULL
  for (i in p_endogen) { # for each lag p_endogen
    for (j in p_exogen) { # for each lag p_exogen
      for (k in s) {
        for (rank in r) {
         for (d_beta in grid_beta) {
          for (d_sigma in grid_sigma) {
          
          pos <- NULL
          model_i <- model
          if (use_discount) {
            model_i[["delta_beta"]] <- d_beta
            model_i[["delta_sigma"]] <- d_sigma
          }
          
          # Build the matrix of regressors from object x
          if (i > 1) {
            pos <- c(pos, 1:(n_endogen * (i - 1)))
          }
          if (i >= 1) {
            model_i[["p"]] <- as.integer(i)
            model_i[["p_endogen"]] <- model_i[["p"]]
          }  
          
          pos <- c(pos, n_endogen * (p_endogen_max - 1) + seq_len(n_exogen * (j - 1)))
          model_i[["p_exogen"]] <- as.integer(j)
          
          if (use_global) {
            pos <- c(pos, n_endogen * (p_endogen_max - 1) + n_exogen * d_exogen_max + seq_len(n_global * (k - 1)))
            model_i[["m_global"]] <- as.integer(n_global)
            model_i[["s_global"]] <- as.integer(k)
          }
          
          # 'm' counts the columns of differences, which is one block less than
          # the lag orders in levels.
          model_i[["m"]] <- model_i[["k_exogen"]] * (model_i[["p_exogen"]] - 1L) +
            model_i[["m_global"]] * max(model_i[["s_global"]] - 1L, 0L)
          model_i[["s"]] <- 1L
          
          if (use_det_ur) {
            pos <- c(pos, n_endogen * (p_endogen_max - 1) + n_exogen * d_exogen_max + n_global * d_global_max + 1:length(det_name_ur))
          }
          
          model_i[["k_beta"]] <- as.integer(n_ect)
          model_i[["rank"]] <- as.integer(rank)
          
          x_i <- NULL
          z <- NULL
          if (length(pos) > 0) {
            # Create data input matrix of the respective model
            x_i <- stats::ts(as.matrix(x[, pos]), class = c("mts", "ts", "matrix"))
            stats::tsp(x_i) <- stats::tsp(temp)
            dimnames(x_i)[[2]] <- x_names[pos]
            z <- kronecker(x_i, diag(1, n_endogen))
          }
          
          if (rank > 0) {
            z <- cbind(matrix(NA_real_, tt * n_endogen, rank * n_endogen) , z)
          }
          
          # If specified add structural data to SUR form
          if (!is.null(y_A0)) {
            z <- cbind(z, y_A0)
          }
          dimnames(z) <- NULL

          # Not kept for a discounted sub-model, which reads the compact
          # regressors: BayesTS refuses the SUR matrix, and it is nine tenths
          # of the sub-model's size.
          if (use_discount) {
            z <- NULL
          }
          
          # Create individual model
          result_i <- list("model" = model_i,
                           "data" = list("original" = list("endogen" = endogen,
                                                           "exogen" = exogen,
                                                           "deterministic" = det_data),
                                         "train" = list("y" = y,
                                                        "w" = ect,
                                                        "x" = x_i,
                                                        "z" = z)))
          
          # Update class of individual model
          class(result_i) <- c("vecxsubmodel", "bvecmodel", "list")
          
          result <- c(result, list(result_i))
          }
         }
        }
      }
    }
  }
  
  
  class(result) <- c("modellist", "list")

  return(result)
}

# The thinning interval a sub-model's sampler keeps one draw in, validated as
# bvartools validates it for its own models: NULL for one, which is how a model
# that keeps every draw says so, and an integer otherwise.
.check_sampler_thin <- function(thin) {
  if (!is.numeric(thin) || length(thin) != 1 || is.na(thin) ||
      thin < 1 || thin != round(thin)) {
    stop("Argument 'thin' must be a single positive integer.")
  }
  if (thin == 1) {
    return(NULL)
  }
  as.integer(thin)
}