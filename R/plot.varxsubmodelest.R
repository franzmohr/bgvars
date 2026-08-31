#' Plotting Draws of a VARX Submodel of a GVAR Model
#' 
#' A plot function for objects of class \code{"varxsubmodelest"} for visual inspection
#' of posterior draws.
#' 
#' @param x an object of class \code{"varxsubmodelest"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param variables character vector of variables that should be plotted. Default is \code{"all"}.
#' Other options are \code{"domestic"}, \code{"foreign"}, \code{"global"}, \code{"deterministic"}
#' and \code{"sigma"}.
#' @param ctry character (optional). Name of the country, which will be shown in the title of the output.
#' @param ... further graphical parameters.
#' 
#' @export
plot.varxsubmodelest <- function(x, ci = 0.95, type = "hist", variables = "all", ctry = NULL, ...) {
  
  orig_par <- par(mar = par("mar"))
  on.exit(par(orig_par))
  
  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }
  
  if (!variables %in% c("all", "domestic", "foreign", "global", "deterministic", "sigma")) {
    stop("Invalid specification of argument 'variables'.")
  }
  
  k_domestic <- x[["model"]][["k_domestic"]]
  k_foreign <- x[["model"]][["k_foreign"]]
  m <- x[["model"]][["m"]]
  n <- x[["model"]][["n"]]
  p_domestic <- x[["model"]][["p_domestic"]]
  p_foreign <- x[["model"]][["p_foreign"]]
  s <- x[["model"]][["s"]]
  tvp <- x[["model"]][["tvp"]]
  
  n_a_domestic <- k_domestic * k_domestic * p_domestic
  n_a_foreign <- k_domestic * k_foreign * (p_foreign + 1)
  n_b <- k_domestic * m * (s + 1)
  
  tt <- NROW(x[["data"]][["y"]])
  tsp_info <- stats::tsp(x[["data"]][["y"]])
  structural <- x[["model"]][["structural"]]
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- dimnames(x[["data"]][["domestic"]])[[2]]
  x_names <- .get_regressor_names_var(x)
  lab_size <- .05
  mar_orig <- graphics::par("mar")
  
  # Title
  title_text <- "Bayesian "
  tvp <- x[["model"]][["tvp"]]
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  sv <- x[["model"]][["error"]] %in% c("sv", "sv+covar")
  if (sv) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VARX submodel")
  
  # Domestic variables ----
  if (p_domestic > 0) {
    
    pos_domestic <- 1:(k_domestic * p_domestic)
    
    .make_new_plot_view(y_names = y_names,
                        x_names = x_names[pos_domestic],
                        title_text = title_text,
                        lab_size = lab_size)
    
    for (i in 1:p_domestic) {
      var_pos <- ((i - 1) * k_domestic * k_domestic) + 1:(k_domestic * k_domestic)
      if (tvp) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["posteriors"]][["a"]][["coeffs"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["a"]][["coeffs"]][, j] == x[["posteriors"]][["a"]][["coeffs"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a"]][["coeffs"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["a"]][["coeffs"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["a"]][["coeffs"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["a"]][["coeffs"]][, j])
            }
          } 
        }
      }
    }
  }
  
  # Foreign variables ----
  pos_foreign <- k_domestic * p_domestic + 1:(k_foreign * (p_foreign + 1))
  
  .make_new_plot_view(y_names = y_names,
                      x_names = x_names[pos_foreign],
                      title_text = title_text,
                      lab_size = lab_size)
  
  for (i in 1:(p_foreign + 1)) {
    var_pos <- n_a_domestic + ((i - 1) * k_domestic * k_foreign) + 1:(k_domestic * k_foreign)
    if (tvp) {
      for (j in var_pos) {
        draws <- .tvpribbon(x[["posteriors"]][["a"]][["coeffs"]], j, ci_low, ci_high)
        stats::tsp(draws) <- tsp_info
        stats::ts.plot(draws, xlab = "")
      }
    } else {
      for (j in var_pos) {
        if (all(x[["posteriors"]][["a"]][["coeffs"]][, j] == x[["posteriors"]][["a"]][["coeffs"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a"]][["coeffs"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["posteriors"]][["a"]][["coeffs"]][, j], plot = TRUE, main = NA)  
          }
          if (type == "trace") {
            stats::ts.plot(x[["posteriors"]][["a"]][["coeffs"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["posteriors"]][["a"]][["coeffs"]][, j])
          }
        }
      }
    }
  }
  
  if (m > 0) {
    stop("implement global")
    
    pos_global <- k_domestic * p_domestic + k_foreign * (p_foreign + 1) + 1:(m * (s + 1))
    
    .make_new_plot_view(y_names = y_names,
                        x_names = x_names[pos_global],
                        title_text = title_text,
                        lab_size = lab_size)
    
    for (i in 1:s) {
      var_pos <- n_a_domestic + n_a_foreign + ((i - 1) * k_domestic * m) + 1:(k_domestic * m)
      if (tvp) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["posteriors"]][["gamma_global"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["gamma_global"]][, j] == x[["posteriors"]][["gamma_global"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["gamma_global"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["gamma_global"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["gamma_global"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["gamma_global"]][, j])
            }
          }
        }
      }
    }
  }
  
  if (n > 0) {
    
    pos_det <- k_domestic * p_domestic + k_foreign * (p_foreign + 1) + m * (s + 1) + 1:n
    
    .make_new_plot_view(y_names = y_names,
                        x_names = x_names[pos_det],
                        title_text = title_text,
                        lab_size = lab_size)
    
    var_pos <- n_a_domestic + n_a_foreign + n_b + 1:(k_domestic * n)
    
    if (tvp) {
      for (j in var_pos) {
        draws <- .tvpribbon(x[["posteriors"]][["a"]][["coeffs"]], j, ci_low, ci_high)
        stats::tsp(draws) <- tsp_info
        stats::ts.plot(draws, xlab = "")
      }
    } else {
      for (j in var_pos) {
        if (all(x[["posteriors"]][["a"]][["coeffs"]][, j] == x[["posteriors"]][["a"]][["coeffs"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a"]][["coeffs"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["posteriors"]][["a"]][["coeffs"]][, j], plot = TRUE, main = NA)  
          }
          if (type == "trace") {
            stats::ts.plot(x[["posteriors"]][["a"]][["coeffs"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["posteriors"]][["a"]][["coeffs"]][, j])
          }
        }
      }
    }
  }
  
  ## Structural ----
  
  if (structural) {
    
    stop("implement structural")
    
    pos_a0 <- k_ect + k_domestic * (p_domestic - 1) + k_foreign * p_foreign + m * s + 1:n_unrestricted
    
    .make_new_plot_view(y_names = y_names,
                        x_names = x_names[pos_det],
                        title_text = title_text,
                        lab_size = lab_size)
    
    var_pos <- n_alpha + n_gamma_domestic + n_gamma_foreign + n_upsilon + 1:(k_domestic * n_unrestricted)
    
    if (tvp) {
      for (j in 1:NCOL(x[["posteriors"]][["a"]][["coeffs"]][[1]])) {
        draws <- .tvpribbon(x[["posteriors"]][["a"]][["coeffs"]], j, ci_low, ci_high)
        if (all(draws[, 1] == draws[1, 1])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
        } else {
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      }
    } else {
      for (j in 1:NCOL(x[["posteriors"]][["a"]][["coeffs"]])) {
        if (all(x[["posteriors"]][["a"]][["coeffs"]][, j] == x[["posteriors"]][["a"]][["coeffs"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a"]][["coeffs"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["posteriors"]][["a"]][["coeffs"]][, j], plot = TRUE, main = NA)   
          }
          if (type == "trace") {
            stats::ts.plot(x[["posteriors"]][["a"]][["coeffs"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["posteriors"]][["a"]][["coeffs"]][, j])
          }
        }
      }
    }
  }
  
  # Sigma ----
  .make_new_plot_view(y_names = y_names,
                      x_names = y_names,
                      title_text = title_text,
                      lab_size = lab_size)
  
  var_pos <- 1:(k_domestic * k_domestic)
  if (sv) {
    for (j in var_pos) {
      draws <- .tvpribbon(x[["posteriors"]][["sigma"]][["coeffs"]], j, ci_low, ci_high)
      if (all(draws[, 1] == draws[1, 1])) {
        graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
      } else {
        stats::tsp(draws) <- tsp_info
        stats::ts.plot(draws, xlab = "") 
      }
    }
  } else {
    for (j in var_pos) {
      if (all(x[["posteriors"]][["sigma"]][["coeffs"]][, j] == x[["posteriors"]][["sigma"]][["coeffs"]][1, j])) {
        graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["sigma"]][["coeffs"]][1, j], adj = 0.5)
      } else {
        if (type == "hist") {
          graphics::hist(x[["posteriors"]][["sigma"]][["coeffs"]][, j], plot = TRUE, main = NA)  
        }
        if (type == "trace") {
          stats::ts.plot(x[["posteriors"]][["sigma"]][["coeffs"]][, j], xlab = "")
        }
        if (type == "boxplot") {
          graphics::boxplot(x[["posteriors"]][["sigma"]][["coeffs"]][, j])
        }
      }
    }
  }
}


