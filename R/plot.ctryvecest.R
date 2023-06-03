#' Plotting Draws of a VECX Submodel of a GVAR Model
#' 
#' A plot function for objects of class \code{"ctryvecest"} for visual inspection
#' of posterior draws.
#' 
#' @param x an object of class \code{"ctryvecest"}, usually, a result of a call to \code{\link{draw_posterior}}.
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
plot.ctryvecest <- function(x, ci = 0.95, type = "hist", variables = "all", ctry = NULL, ...) {
  
  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }
  
  if (!variables %in% c("all", "domestic", "foreign", "deterministic", "sigma")) {
    stop("Invalid specification of argument 'variables'.")
  }
  
  x <- .create_pi_matrices(x)
  
  names_domestic <- x[["model"]][["domestic"]][["variables"]]
  k_domestic <- length(x[["model"]][["domestic"]][["variables"]])
  p_domestic <- x[["model"]][["domestic"]][["lags"]]
  names_foreign <- x[["model"]][["foreign"]][["variables"]]
  k_foreign <- length(x[["model"]][["foreign"]][["variables"]])
  p_foreign <- x[["model"]][["foreign"]][["lags"]]
  global <- !is.null(x[["model"]][["global"]][["variables"]])
  if (global) {
    names_global <- x[["model"]][["global"]][["variables"]]
    k_global <- length(x[["model"]][["global"]][["variables"]])
    s_global <- x[["model"]][["global"]][["lags"]] 
  } else {
    k_global <- 0
    s_global <- 0
  }
  names_det_r <- x[["model"]][["deterministic"]][["restricted"]]
  k_det_r <- length(names_det_r)
  names_det_ur <- x[["model"]][["deterministic"]][["unrestricted"]]
  k_det_ur <- length(names_det_ur)
  r <- x[["model"]][["rank"]]
  
  tt <- NROW(x[["data"]][["Y"]])
  tsp_info <- stats::tsp(x[["data"]][["Y"]])
  structural <- x[["model"]][["structural"]]
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- paste0("d.", dimnames(x[["data"]][["domestic"]])[[2]])
  x_names <- .get_regressor_names_vec(x)
  lab_size <- .05
  mar_orig <- graphics::par("mar")
  
  # Calculate number of columns in final layout
  n_tot <- 0
  if ("pi_domestic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_domestic
  }
  if ("pi_foreign" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_foreign
  }
  if ("pi_global" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_global
  }
  if ("pi_deterministic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + length(names_det_r)
  }
  if ("gamma_domestic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_domestic * (p_domestic - 1)
  }
  if ("gamma_foreign" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_foreign * p_foreign
  }
  if ("gamma_global" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_global * s_global
  }
  if ("gamma_deterministic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + length(names_det_ur)
  }
  if (!is.null(x[["posteriors"]][["a0"]])) {
    n_tot <- n_tot + k_domestic
  }
  if ("sigma" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_domestic
  }
  
  if (variables == "domestic") {
    n_tot <- k_domestic * p_domestic
    if (p_domestic > 1) {
      x_names <- x_names[c(1:k_domestic, k_domestic + k_foreign + k_global + k_det_r + 1:(k_domestic * (p_domestic - 1)))] 
    } else {
      x_names <- x_names[1:k_domestic] 
    }
  }
  if (variables == "foreign") {
    n_tot <- k_foreign * (p_foreign + 1)
    x_names <- x_names[c(k_domestic + 1:k_foreign, k_domestic * p_domestic + k_foreign + k_global + k_det_r + 1:(k_foreign * (p_foreign)))]
  }
  if (variables == "global") {
    n_tot <- k_global * (s_global + 1)
    x_names <- x_names[c(k_domestic + k_foreign + 1:k_global,
                         k_domestic * p_domestic + k_foreign * (p_foreign + 1) + k_global + k_det_r + 1:(k_global * (s_global)))]
  }
  if (variables == "deterministic") {
    n_tot <- k_det_r + k_det_ur
    pos_temp <- NULL
    if (k_det_r > 0) {
      pos_temp <- k_domestic + k_foreign + k_global + 1:k_det_r
    }
    if (k_det_ur > 0) {
      pos_temp <- c(pos_temp, k_domestic * p_domestic + k_foreign * (p_foreign + 1) + k_global * (s_global + 1) + k_det_r + 1:k_det_ur)
    }
    x_names <- x_names[pos_temp]
  }
  if (variables == "sigma") {
    n_tot <- k_domestic
  }
  
  # Define layout
  mat <- matrix(NA_integer_,
                k_domestic + 2 , # k_domestic + title row + column names
                n_tot + 1) # n_tot + row names
  mat[1, ] <- 1
  mat[-1, 1] <- c(0, 2:(k_domestic + 1))
  mat[2, -1] <- (k_domestic + 1) + 1:n_tot
  mat[-(1:2), -1] <- matrix(1:(k_domestic * n_tot) + k_domestic + n_tot + 1, k_domestic, n_tot)
  graphics::layout(mat,
                   widths = c(lab_size, rep((1 - lab_size) / n_tot, n_tot)),
                   heights = c(.07, lab_size, rep((1 - lab_size) / k_domestic, k_domestic)))
  
  # Title
  title_text <- "Bayesian "
  tvp <- x[["model"]][["tvp"]]
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (x[["model"]][["sv"]]) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VECX submodel")
  if (!is.null(ctry)) {
    title_text <- paste0(title_text, " - ", ctry)
  }
  
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new(); graphics::text(0.5, 0.5, labels = title_text, cex = 1.5)
  # Fill rows
  graphics::par(mar = c(3, 0, 0, 0))
  for (j in y_names) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  # Fill columns
  graphics::par(mar = c(0, 0, 0, 0))
  if (variables %in% c("all", "domestic", "foreign", "global", "deterministic")) {
    for (j in x_names) {
      graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
    } 
  }
  if (variables %in% c("all", "sigma")) {
    for (j in y_names) {
      graphics::plot.new(); graphics::text(0.5, 0.5, labels = paste0("Sigma\n", j), adj = 0.5)
    }  
  }
  
  graphics::par(mar = c(3, 2.1, .5, 1))
  
  if (variables %in% c("all", "domestic")) {
    if (!is.null(x[["posteriors"]][["pi_domestic"]])) {
      var_pos <- 1:(k_domestic * k_domestic)
      if (tvp) {
        for (j in var_pos) {
          # draws <- .tvpribbon(x[["A"]], j, ci_low, ci_high)
          # stats::tsp(draws) <- tsp_info
          # stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["pi_domestic"]][, j] == x[["posteriors"]][["pi_domestic"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["pi_domestic"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["pi_domestic"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["pi_domestic"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["pi_domestic"]][, j])
            }
          } 
        }
      }
    }    
  }
  
  if (variables %in% c("all", "foreign")) {
    if (!is.null(x[["posteriors"]][["pi_foreign"]])) {
      var_pos <- 1:(k_domestic * k_foreign)
      if (tvp) {
        for (j in var_pos) {
          # draws <- .tvpribbon(x[["A"]], j, ci_low, ci_high)
          # stats::tsp(draws) <- tsp_info
          # stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["pi_foreign"]][, j] == x[["posteriors"]][["pi_foreign"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["pi_foreign"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["pi_foreign"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["pi_foreign"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["pi_foreign"]][, j])
            }
          } 
        }
      }
    }
  }
  
  if (variables %in% c("all", "global")) {
    if (!is.null(x[["posteriors"]][["pi_global"]])) {
      var_pos <- 1:(k_domestic * k_global)
      if (tvp) {
        for (j in var_pos) {
          # draws <- .tvpribbon(x[["A"]], j, ci_low, ci_high)
          # stats::tsp(draws) <- tsp_info
          # stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["pi_global"]][, j] == x[["posteriors"]][["pi_global"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["pi_global"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["pi_global"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["pi_global"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["pi_global"]][, j])
            }
          } 
        }
      }
    }
  }
  
  if (variables %in% c("all", "deterministic")) {
    if (!is.null(x[["posteriors"]][["pi_deterministic"]])) {
      var_pos <- 1:(k_domestic * length(names_det_r))
      if (tvp) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["pi_deterministic"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["pi_deterministic"]][, j] == x[["posteriors"]][["pi_deterministic"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["pi_deterministic"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["pi_deterministic"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["pi_deterministic"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["pi_deterministic"]][, j])
            }
          } 
        }
      }
    }
  }
  
  if (variables %in% c("all", "domestic")) {
    if (!is.null(x[["posteriors"]][["gamma_domestic"]])) {
      if (p_domestic > 1) {
        for (i in 1:(p_domestic - 1)) {
          var_pos <- ((i - 1) * k_domestic * k_domestic) + 1:(k_domestic * k_domestic)
          if (tvp) {
            for (j in var_pos) {
              draws <- .tvpribbon(x[["gamma_domestic"]], j, ci_low, ci_high)
              stats::tsp(draws) <- tsp_info
              stats::ts.plot(draws, xlab = "")
            }
          } else {
            for (j in var_pos) {
              if (all(x[["posteriors"]][["gamma_domestic"]][, j] == x[["posteriors"]][["gamma_domestic"]][1, j])) {
                graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["gamma_domestic"]][1, j], adj = 0.5)
              } else {
                if (type == "hist") {
                  graphics::hist(x[["posteriors"]][["gamma_domestic"]][, j], plot = TRUE, main = NA)  
                }
                if (type == "trace") {
                  stats::ts.plot(x[["posteriors"]][["gamma_domestic"]][, j], xlab = "")
                }
                if (type == "boxplot") {
                  graphics::boxplot(x[["posteriors"]][["gamma_domestic"]][, j])
                }
              } 
            }
          }
        } 
      }
    }
  }
  
  if (variables %in% c("all", "foreign")) {
    if (!is.null(x[["posteriors"]][["gamma_foreign"]])) {
      for (i in 1:p_foreign) {
        var_pos <- ((i - 1) * k_domestic * k_foreign) + 1:(k_domestic * k_foreign)
        if (tvp) {
          for (j in var_pos) {
            draws <- .tvpribbon(x[["gamma_foreign"]], j, ci_low, ci_high)
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        } else {
          for (j in var_pos) {
            if (all(x[["posteriors"]][["gamma_foreign"]][, j] == x[["posteriors"]][["gamma_foreign"]][1, j])) {
              graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["gamma_foreign"]][1, j], adj = 0.5)
            } else {
              if (type == "hist") {
                graphics::hist(x[["posteriors"]][["gamma_foreign"]][, j], plot = TRUE, main = NA)  
              }
              if (type == "trace") {
                stats::ts.plot(x[["posteriors"]][["gamma_foreign"]][, j], xlab = "")
              }
              if (type == "boxplot") {
                graphics::boxplot(x[["posteriors"]][["gamma_foreign"]][, j])
              }
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "global")) {
    if (!is.null(x[["posteriors"]][["gamma_global"]])) {
      for (i in 1:s_global) {
        var_pos <- ((i - 1) * k_domestic * k_global) + 1:(k_domestic * k_global)
        if (tvp) {
          for (j in var_pos) {
            draws <- .tvpribbon(x[["gamma_global"]], j, ci_low, ci_high)
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
  }
  
  if (variables %in% c("all", "deterministic")) {
    if (!is.null(x[["posteriors"]][["gamma_deterministic"]])) {
      if (tvp) {
        for (j in 1:NCOL(x[["posteriors"]][["gamma_deterministic"]][[1]])) {
          draws <- .tvpribbon(x[["posteriors"]][["gamma_deterministic"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in 1:NCOL(x[["posteriors"]][["gamma_deterministic"]])) {
          if (all(x[["posteriors"]][["gamma_deterministic"]][, j] == x[["posteriors"]][["gamma_deterministic"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["gamma_deterministic"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["gamma_deterministic"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["gamma_deterministic"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["gamma_deterministic"]][, j])
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all")) {
    if (!is.null(x[["posteriors"]][["a0"]])) {
      if (tvp) {
        for (j in 1:NCOL(x[["posteriors"]][["a0"]][[1]])) {
          draws <- .tvpribbon(x[["posteriors"]][["a0"]], j, ci_low, ci_high)
          if (all(draws[, 1] == draws[1, 1])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
          } else {
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        }
      } else {
        for (j in 1:NCOL(x[["posteriors"]][["a0"]])) {
          if (all(x[["posteriors"]][["a0"]][, j] == x[["posteriors"]][["a0"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a0"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["a0"]][, j], plot = TRUE, main = NA)   
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["a0"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["a0"]][, j])
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "sigma")) {
    if (!is.null(x[["posteriors"]][["sigma"]])) {
      var_pos <- 1:(k_domestic * k_domestic)
      if (x[["model"]][["structural"]]) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["posteriors"]][["sigma"]], j, ci_low, ci_high)
          if (all(draws[, 1] == draws[1, 1])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
          } else {
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "") 
          }
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["sigma"]][, j] == x[["posteriors"]][["sigma"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["sigma"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["sigma"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["sigma"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["sigma"]][, j])
            }
          }
        }
      }
    }
  }
  
  graphics::par(mar = mar_orig)
  graphics::layout(matrix(1))
}


