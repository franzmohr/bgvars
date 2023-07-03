#' Test Statistics for the Selection of Submodels of GVAR Model
#'
#' Calculates test statistics for the the selection of country-specific VARX or 
#' VECX models of a GVAR model.
#'
#' @param object an object of class \code{"bgvarest"} or \code{"bgvecest"}, usually,
#' a result of a call to \code{\link{draw_posterior.gvarsubmodels}} or 
#' \code{\link{draw_posterior.gvecsubmodels}}, respectively.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The log-likelihood for the calculation of the information criteria is obtained by
#' \deqn{LL = \frac{1}{R} \sum_{i = 1}^{R} \left( \sum_{t = 1}^{T} -\frac{K^{dom}}{2} \ln 2\pi - \frac{1}{2} \ln |\Sigma_t^{(i)}| -\frac{1}{2} (u_t^{{(i)}\prime} (\Sigma_t^{(i)})^{-1} u_t^{(i)} \right)},
#' where \eqn{u_t = y_t - \mu_t}.
#' 
#' For VAR models the Akaike, Bayesian and Hannan–Quinn (HQ) information criteria are calculated as
#' \deqn{AIC = 2 (K^{d}p^{d} + K^{f}p^{f} + Ms + N) - 2 LL},
#' \deqn{BIC = (K^{d}p^{d} + K^{f}p^{f} + Ms + N) ln(T) - 2 LL} and 
#' \deqn{HQ = 2  (K^{d}p^{d} + K^{f}p^{f} + Ms + N) ln(ln(T)) - 2 LL}, respectively,
#' where \eqn{K^{d}} is the number of endogenous domestic variables, \eqn{p^{d}} the number of lags of endogenous domestic variables,
#' \eqn{K^{f}} is the number of foreign variables, \eqn{p^{f}} the number of lags of foreign variables,
#' \eqn{M} the number of global variables, \eqn{s} the number of lags of global variables,
#' \eqn{N} the number of deterministic terms and \eqn{T} the number of observations.
#'
#' @return A list.
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
#' # Obtain test statistics
#' tests <- submodel_test_statistics(object, ic = "BIC", select = "order")
#' 
#'
#' @export
submodel_test_statistics <- function(object, ...){
  
  n_models <- length(object)
  
  teststats <- data.frame(ctry = rep(NA, n_models),
                          p_domestic = rep(NA, n_models),
                          p_foreign = rep(NA, n_models),
                          s = rep(NA, n_models),
                          r = rep(NA, n_models),
                          LL = rep(NA, n_models),
                          AIC = rep(NA, n_models),
                          BIC = rep(NA, n_models),
                          HQ = rep(NA, n_models),
                          stringsAsFactors = FALSE)
  
  loglik <- list()
  
  for (i in 1:n_models) {
    
    teststats[i, "ctry"] <- names(object)[i]
    
    # Skip tests if posterior simulation was not successful
    if (!is.null(object[[i]][["error"]])) {
      if (object[[i]][["error"]]) {
        next 
      }
    }
    
    structural <- object[[i]][["model"]][["structural"]]
    tvp <- object[[i]][["model"]][["tvp"]]
    sv <- object[[i]][["model"]][["sv"]]
    tt <- NROW(object[[i]][["data"]][["Y"]])
    k_domestic <- length(object[[i]][["model"]][["domestic"]][["variables"]])
    p_domestic <- object[[i]][["model"]][["domestic"]][["lags"]]
    k_foreign <- length(object[[i]][["model"]][["foreign"]][["variables"]])
    p_foreign <- object[[i]][["model"]][["foreign"]][["lags"]]
    global <- !is.null(object[[i]][["model"]][["global"]])
    if (global) {
      k_global <- length(object[[i]][["model"]][["global"]][["variables"]])
      s <- object[[i]][["model"]][["global"]][["lags"]]
    }
    
    teststats[i, "p_domestic"] <- p_domestic
    teststats[i, "p_foreign"] <- p_foreign 
    if (global) {
      teststats[i, "s"] <- s 
    }
    if (!is.null(object[[i]][["model"]][["rank"]])) {
      teststats[i, "r"] <- object[[i]][["model"]][["rank"]]
    }
    
    if (tvp) {
      temp_pars <- list()
      length(temp_pars) <- tt
    } else {
      temp_pars <- NULL
    }
    x <- NULL
    
    if ("ctryvarest" %in% class(object[[i]])) {
      
      x <- t(object[[i]][["data"]][["Z"]])
      tot_pars <- NCOL(object[[i]][["data"]][["Z"]])
      
      vars <- c("domestic", "foreign", "global", "deterministic")
      for (j in vars) {
        if (!is.null(object[[i]][["posteriors"]][[j]])) {
          if (is.list(object[[i]][["posteriors"]][[j]])) {
            for (period in 1:tt) {
              temp_pars[[period]] <- cbind(temp_pars[[period]], object[[i]][["posteriors"]][[j]][[period]]) 
            }
          } else {
            temp_pars <- cbind(temp_pars, object[[i]][["posteriors"]][[j]]) 
          }
        }
      }
    }
    
    if ("ctryvecest" %in% class(object[[i]])) {
      
      object[[i]] <- .create_pi_matrices(object[[i]])
      
      x <- t(object[[i]][["data"]][["X"]])
      tot_pars <- object[[i]][["model"]][["rank"]] + nrow(x)
      if (object[[i]][["model"]][["rank"]] > 0) {
        x <- rbind(t(object[[i]][["data"]][["W"]]), x)
      }
      
      vars <- c("pi_domestic", "pi_foreign", "pi_global", "pi_deterministic",
                "gamma_domestic", "gamma_foreign", "gamma_global", "gamma_deterministic")
      for (j in vars) {
        if (!is.null(object[[i]][["posteriors"]][[j]])) {
          if (is.list(object[[i]][["posteriors"]][[j]])) {
            for (period in 1:tt) {
              temp_pars[[period]] <- cbind(temp_pars[[period]], object[[i]][["posteriors"]][[j]][[period]]) 
            }
          } else {
            temp_pars <- cbind(temp_pars, object[[i]][["posteriors"]][[j]]) 
          }
        }
      }
      
    }
    
    if (tvp) {
      draws <- nrow(temp_pars[[1]])
    } else {
      draws <- nrow(temp_pars) 
    }
    loglik[[i]] <- matrix(NA, tt, draws)
    y <- t(object[[i]][["data"]][["Y"]])
    u <- y * 0
    if (sv) {
      sigma <- matrix(NA_real_, k_domestic * tt, k_domestic)
    } else {
      sigma <- matrix(NA_real_, k_domestic, k_domestic)
    }
    
    for (j in 1:draws) {
      # Residuals
      if (tvp) {
        for (period in 1:tt) {
          if (structural) {
            A0 <- matrix(object[[i]][["posteriors"]][["a0"]][[period]][j, ], k_domestic)
          } else {
            A0 <- diag(k_domestic)
          }
          u[, period] <- y[, period] - matrix(temp_pars[[period]][j, ], k_domestic) %*% x[, period] 
        }
      } else {
        if (structural) {
          A0 <- matrix(object[[i]][["posteriors"]][["a0"]][j, ], k_domestic)
        } else {
          A0 <- diag(k_domestic)
        }
        u <- A0 %*% y - matrix(temp_pars[j, ], k_domestic) %*% x 
      }
      
      if (sv) {
        for (period in 1:tt) {
          sigma[(period - 1) * k_domestic + 1:k_domestic,] <- matrix(object[[i]][["posteriors"]][["sigma"]][[period]][j,], k_domestic)
        }
      } else {
        sigma <- matrix(object[[i]][["posteriors"]][["sigma"]][j,], k_domestic)
      }
      
      # Calculate log-likelihood
      loglik[[i]][, j] <- bvartools::loglik_normal(u, sigma)
    }
    
    ll <- sum(rowMeans(loglik[[i]]))
    teststats[i, "LL"] <- ll
    teststats[i, "AIC"] <- 2 * tot_pars - 2 * ll
    teststats[i, "BIC"] <- log(tt) * tot_pars - 2 * ll
    teststats[i, "HQ"] <- 2 * log(log(tt)) * tot_pars - 2 * ll
  }
  
  # Omit unnecessary columns
  teststats <- teststats[, which(!apply(teststats, 2, function(x) {all(is.na(x))}))]
  
  # Final output with one list element per country
  result <- NULL
  ctry_names <- unique(teststats[, "ctry"])
  for (i in 1:length(ctry_names)) {
    result[[i]] <- teststats[teststats[, "ctry"] == ctry_names[i], -1]
    rownames(result[[i]]) <- NULL
    names(result)[i] <- ctry_names[i]
  }
  names(loglik) <- teststats[, "ctry"]
  
  result <- list(teststats = teststats,
                 loglik = loglik)
  
  return(result)
}