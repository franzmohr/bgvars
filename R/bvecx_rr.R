.bvecx_rr <- function(object, iterations, burnin, thin){
  #### Specifications ####
  #rm(list = ls()[-which(ls() %in% c("object", "iterations", "burnin", "thin"))])
  estimation_data <- .gen_vecx(object)
  y <- estimation_data$y
  ect <- estimation_data$ect
  x <- estimation_data$x
  ect_x <- rbind(ect, x)
  tt <- dim(y)[2]
  
  # Specify rank
  r <- object$model$cointegration$rank
  use_rr <- r > 0
  
  # Number of variables and lags
  k_domestic <- estimation_data$domestic$dim
  p_domestic <- estimation_data$domestic$lag
  k_foreign <- estimation_data$foreign$dim
  p_foreign <- estimation_data$foreign$lag
  
  k_ect <- nrow(ect)
  k_x <- nrow(x)
  n_domestic <- k_domestic * k_domestic
  n_foreign <- k_domestic * k_foreign 
  n_ect <- k_domestic * k_ect
  n_x <- k_domestic * k_x
  n_alpha <- k_domestic * r
  n_beta <- k_ect * r
  
  global <- !is.null(estimation_data$global)
  k_global <- 0
  n_global <- 0
  if (global) {
    k_global <- estimation_data$global$dim
    p_global <- estimation_data$global$lag
    n_global <- k_domestic * k_global
  }
  
  k_res <- estimation_data$deterministics$restricted
  n_res <- k_domestic * k_res
  k_unres <- estimation_data$deterministics$unresticted
  n_unres <- k_domestic * k_unres
  
  #### Priors ####
  # Non-cointegration
  gamma_mu_prior <- object$priors$coefficients$mu
  gamma_v_i_prior <- object$priors$coefficients$v_i
  
  if (use_rr) {
    v_i <- object$priors$cointegration$v_i
    p_tau_i <- object$priors$cointegration$p_tau_i 
  }
  
  # Errors
  sigma_df_prior <- object$priors$sigma$df
  sigma_df_post <- sigma_df_prior + tt
  sigma_scale_prior <- object$priors$sigma$scale
  
  #### Variable selection ####
  var_select <- object$model$variable_selection$type != "none"
  gamma_update <- FALSE
  if (var_select) {
    var_select_type <- object$model$variable_selection$type
    gamma_select_draws <- iterations - burnin
    gamma_update <- TRUE # Is TRUE during sampling unless 'draws' was specified
    gamma_lambda <- diag(1, n_x) # Initial inclusion matrix
    gamma_prob_prior <- object$priors$variable_selection$inprior # Prior inclusion probabilities
    gamma_include <- object$model$variable_selection$include # Indicator for which coeffs VS should be done
    if (var_select_type == "bvs") {
      sur <- TRUE
    }
    if (var_select_type == "ssvs") {
      tau0 <- matrix(object$priors$variable_selection$tau0[n_ect + 1:n_x,])
      tau1 <- matrix(object$priors$variable_selection$tau1[n_ect + 1:n_x,])
      gamma_v_i_prior <- solve(diag(c(tau1^2)))
    }
    if (!is.null(object$model$variable_selection$threshol)) {
      vs_threshold <- object$model$variable_selection$threshold
      if (!is.null(object$model$variable_selection$draws)) {
        gamma_select_draws <- object$model$variable_selection$draws
      } else {
        gamma_select_draws <- floor((iterations - burnin) / 2)
      }
      if (gamma_select_draws > iterations - burnin) {
        stop("Number of variable selection draws may not be higher than the number of post-burn-in draws.")
      }
    }
  }
  
  #### Initial values ####
  if (use_rr) {
    ols <- tcrossprod(y, ect_x) %*% solve(tcrossprod(ect_x)) 
    u <- y - ols %*% ect_x
  } else {
    ols <- tcrossprod(y, x) %*% solve(tcrossprod(x)) 
    u <- y - ols %*% x
  }
  sigma <- tcrossprod(u) / (tt - nrow(x))
  sigma_i <- solve(sigma)
  rm(list = c("u", "ols"))
  
  y_vec <- matrix(y)
  z <- kronecker(t(x), diag(1, k_domestic))
  
  if (var_select) {
    z_bvs <- z
  }
  
  pos_gamma <- 1:n_x
  
  if (use_rr) {
    beta <- matrix(0, k_ect, r)
    beta[1:r, 1:r] <- diag(1, r)
  }
  
  #### Data container ####
  store <- iterations - burnin
  
  if (use_rr) {
    draws_pi <- matrix(NA, n_ect, store)  
  }
  draws_coef <- matrix(NA, n_x, store) 
  draws_sigma <- matrix(NA, n_domestic, store)
  draws_ll <- matrix(NA, tt, store)
  
  if (var_select) {
    draws_lambda <- matrix(NA, n_x, gamma_select_draws)
  }
  
  #### Gibbs sampler ####
  for (draw in 1:iterations) {
    if (var_select & gamma_update) {
      z <- z_bvs %*% gamma_lambda
    }
    
    if (use_rr) {
      temp <- post_coint_kls_sur(y = y, beta = beta, w = ect, sigma_i = sigma_i,
                                 v_i = v_i, p_tau_i = p_tau_i, g_i = sigma_i,
                                 x = z, gamma_mu_prior = gamma_mu_prior, gamma_v_i_prior = gamma_v_i_prior)
      Pi <- temp$Pi
      gamma <- temp$Gamma
      beta <- temp$beta 
    } else {
      gamma <- post_normal_sur(y = y, z = z, sigma_i = sigma_i,
                               a_prior = gamma_mu_prior, v_i_prior = gamma_v_i_prior)
    }
    
    if (var_select) {
      if (var_select_type == "ssvs" & gamma_update) {
        temp_ssvs <- ssvs(gamma, tau0, tau1, gamma_prob_prior, gamma_include)
        gamma_lambda <- diag(c(temp_ssvs$lambda))
        gamma_v_i_prior <- temp_ssvs$V_i
      }
      
      if (var_select_type == "bvs") {
        if (gamma_update) {
          gamma_lambda <- bvs(y, z_bvs, gamma, gamma_lambda, sigma_i, gamma_prob_prior, gamma_include) 
        }
        gamma <- gamma_lambda %*% gamma
      }
    }
    
    #### Errors ####
    if (use_rr) {
      y_star <- matrix(y_vec - matrix(Pi %*% ect) - z %*% gamma, k_domestic)
    } else {
      y_star <- matrix(y_vec - z %*% gamma, k_domestic)
    }
    
    sigma_scale_i_post <- solve(sigma_scale_prior + tcrossprod(y_star))
    sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_i_post)[,, 1], k_domestic)
    sigma <- solve(sigma_i)
    
    #### Save parameters ####
    if (draw > burnin) {
      pos_draw <- draw - burnin
      
      draws_sigma[, pos_draw] <- sigma
      
      if (use_rr) {
        draws_pi[, pos_draw] <- Pi
      }
      
      draws_coef[, pos_draw] <- gamma
      
      draws_ll[, pos_draw] <- loglik_normal(y_star, sigma)
      
      if (var_select) {
        if (pos_draw <= gamma_select_draws) {
          draws_lambda[, pos_draw] <- diag(gamma_lambda)
          if (pos_draw == gamma_select_draws & draw < iterations) {
            diag(gamma_lambda) <- as.numeric(rowMeans(draws_lambda) >= vs_threshold)
            if (var_select_type == "bvs") {
              z <- z_bvs %*% gamma_lambda 
            }
            if (var_select_type == "ssvs") {
              diag(gamma_v_i_prior)[diag(gamma_lambda) == 0] <- 1 / (tau0[diag(gamma_lambda) == 0,]^2)
              diag(gamma_v_i_prior)[diag(gamma_lambda) == 1] <- 1 / (tau1[diag(gamma_lambda) == 1,]^2)
            }
            gamma_update <- FALSE
          }
        }
      }
      
      #### Print draws ####
      if (FALSE) {
        if (draw %in% floor(seq(from = burnin + 10, to = iterations, length.out = 100))) {
          if (pos_draw >= 10) {
            
            graphics::par(mfcol = c(2, 2))
            temp_pos <- round(.5 * pos_draw):pos_draw
            
            imp <- "d.s.y"
            res <- "d.y"
            pos_var <- k_domestic * (which(dimnames(x)[[1]] == imp) - 1) + which(dimnames(y)[[1]] == res)
            
            stats::plot.ts(draws_coef[pos_var, temp_pos], ylab = "contemporaneous")
            
            
            res <- "d.y"
            if (n_unres > 0) {
              pos_var <- n_x - n_unres + which(dimnames(y)[[1]] == res)
              
              plot.ts(draws_coef[pos_var, temp_pos], ylab = "deterministic")
              
            } else {
              if (use_rr) {
                if (n_res > 0) {
                  pos_var <- n_ect - n_res + which(dimnames(y)[[1]] == res)
                  
                  plot.ts(draws_pi[pos_var, temp_pos], ylab = "deterministic restr")
                  
                } 
              }
            }
            
            imp <- "d.Dp"
            res <- "d.Dp"
            pos_var <- k_domestic * (which(dimnames(y)[[1]] == imp) - 1) + which(dimnames(y)[[1]] == res)
            
            plot.ts(draws_sigma[pos_var, temp_pos], ylab = "covariance") 
            
            
            res <- "d.y"
            pos_var <- k_domestic * (which(dimnames(y)[[1]] == res) - 1) + which(dimnames(y)[[1]] == res)
            
            plot.ts(draws_sigma[pos_var, temp_pos], ylab = "variance")
            
            graphics::par(mfcol = c(1, 1))
          }
        } 
      }
    }
  }
  
  #### Thinning ####
  if (var_select & !is.null(object$model$variable_selection$threshold)) {
    pos_thin <- seq(from = gamma_select_draws + thin, to = store, by = thin)
    burnin <- burnin + gamma_select_draws
  } else {
    pos_thin <- seq(from = thin, to = store, by = thin)
  }
  
  #### Combine coefficients ####
  mc_start <- burnin + thin
  mc_thin <- thin
  
  # Cointegration parameters
  if (use_rr) {
    object$coefficients <- list(pi_domestic = coda::mcmc(t(matrix(draws_pi[1:n_domestic, pos_thin], n_domestic)), start = mc_start, thin = mc_thin),
                                pi_foreign = coda::mcmc(t(matrix(draws_pi[n_domestic + 1:n_foreign, pos_thin], n_foreign)), start = mc_start, thin = mc_thin))
    if (global) {
      object$coefficients$pi_global <- coda::mcmc(t(matrix(draws_pi[n_domestic + n_foreign + 1:n_global, pos_thin], n_global)), start = mc_start, thin = mc_thin)
    } 
    if (k_res > 0) {
      object$coefficients$c_res <- coda::mcmc(t(matrix(draws_pi[n_domestic + n_foreign + n_global + 1:n_res, pos_thin], n_res)), start = mc_start, thin = mc_thin)
    }
  }
  
  # Non-cointegration Parameters
  tot_domestic <- 0
  if (p_domestic > 0) {
    tot_domestic <- n_domestic * p_domestic
    object$coefficients$gamma_domestic <- coda::mcmc(t(matrix(draws_coef[1:tot_domestic, pos_thin], tot_domestic)), start = mc_start, thin = mc_thin) 
  }
  
  tot_foreign <- n_foreign * (p_foreign + 1)
  object$coefficients$gamma_foreign <- coda::mcmc(t(matrix(draws_coef[tot_domestic + 1:tot_foreign, pos_thin], tot_foreign)), start = mc_start, thin = mc_thin)  
  
  tot_global <- 0
  if (global) {
    tot_global <- n_global * (p_global + 1)
    object$coefficients$gamma_global <- coda::mcmc(t(matrix(draws_coef[tot_domestic + tot_foreign + 1:tot_global, pos_thin], tot_global)), start = mc_start, thin = mc_thin) 
  }
  
  if (k_unres > 0) {
    tot_unres <- k_domestic * k_unres
    object$coefficients$c_unres <- coda::mcmc(t(matrix(draws_coef[tot_domestic + tot_foreign + tot_global + 1:tot_unres, pos_thin], tot_unres)), start = mc_start, thin = mc_thin) 
  }
  
  object$coefficients$sigma <- coda::mcmc(t(matrix(draws_sigma[, pos_thin], n_domestic)), start = mc_start, thin = mc_thin) 
  
  if (var_select) {
    lambda_mean <- matrix(rowMeans(draws_lambda), k_domestic)
    dimnames(lambda_mean) <- list(dimnames(y)[[1]], dimnames(x)[[1]])
    object$coefficients$lambda_mean <- lambda_mean
  }
  
  tot_pars <- n_alpha + n_x
  draws_ll <- draws_ll[, pos_thin]
  LL <- sum(log(rowMeans(exp(draws_ll))))
  object$teststats <- c("loglik" = LL,
                        "AIC" = 2 * tot_pars - 2 * LL,
                        "BIC" = log(tt) * tot_pars - 2 * LL,
                        "HQ" = 2 * log(log(tt)) * tot_pars - 2 * LL)
  
  #### VEC to VAR ####
  draws <- nrow(object$coefficients$gamma_foreign) 
  
  
  if (p_domestic > 0) {
    p_domestic <- p_domestic + 1
    W <- diag(-1, k_domestic * p_domestic)
    W[1:k_domestic, 1:k_domestic] <- diag(1, k_domestic)
    W[-(1:k_domestic), -(k_domestic * (p_domestic - 1) + 1:k_domestic)] <- W[-(1:k_domestic),-(k_domestic * (p_domestic - 1) + 1:k_domestic)] + diag(k_domestic * (p_domestic - 1))
    J <- matrix(0, k_domestic, k_domestic * p_domestic)
    J[1:k_domestic, 1:k_domestic] <- diag(1, k_domestic)
    
    temp <- matrix(NA, k_domestic^2 * p_domestic, draws)
    use_rr_temp <- matrix(0, k_domestic, k_domestic)
    for (draw in 1:draws) {
      if (use_rr) {
        use_rr_temp <- matrix(object$coefficients$pi_domestic[draw, ], k_domestic)
      }
      temp[, draw] <- cbind(use_rr_temp,
                            matrix(object$coefficients$gamma_domestic[draw, ], k_domestic)) %*% W + J
    } 
    
  } else {
    p_domestic <- p_domestic + 1
    temp <- matrix(0, n_domestic, draws)
    if (use_rr) {
      for (draw in 1:draws) {
        temp[, draw] <- matrix(matrix(object$coefficients$pi_domestic[draw, ], k_domestic) + diag(1, k_domestic), k_domestic)
      } 
    }
  }
  
  object$coefficients$a_domestic <- coda::mcmc(t(temp), start = mc_start, thin = mc_thin)
  rm(temp) 
  
  
  W <- diag(-1, k_foreign * (p_foreign + 2))
  W[1:k_foreign, 1:k_foreign] <- 0
  W[1:k_foreign, k_foreign + 1:k_foreign] <- diag(1, k_foreign)
  W[-(1:k_foreign), 1:(k_foreign * (p_foreign + 1))] <- W[-(1:k_foreign), 1:(k_foreign * (p_foreign + 1))] + diag(1, k_foreign * (p_foreign + 1))
  
  temp <- matrix(NA, k_domestic * k_foreign * (p_foreign + 2), draws)
  rr_temp <- matrix(0, k_domestic, k_foreign)
  for (draw in 1:draws){
    if (use_rr) {
      rr_temp <- matrix(object$coefficients$pi_foreign[draw, ], k_domestic)
    }
    temp[, draw] <- cbind(rr_temp,
                          matrix(object$coefficients$gamma_foreign[draw, ], k_domestic)) %*% W
  }
  object$coefficients$a_foreign <- coda::mcmc(t(temp), start = mc_start, thin = mc_thin)
  rm(temp) 
  
  if (global) {
    W <- diag(-1, k_global * (p_global + 2))
    W[1:k_global, 1:k_global] <- 0
    W[1:k_global, k_global + 1:k_global] <- diag(1, k_global)
    W[-(1:k_global), 1:(k_global * (p_global + 1))] <- W[-(1:k_global), 1:(k_global * (p_global + 1))] + diag(1, k_global * (p_global + 1))
    
    temp <- matrix(NA, k_domestic * k_global * (p_global + 2), draws)
    rr_temp <- matrix(0, k_domestic, k_global)
    for (draw in 1:draws){
      if (use_rr) {
        use_rr_temp <- matrix(object$coefficients$pi_global[draw, ], k_domestic)
      }
      temp[, draw] <- cbind(rr_temp,
                            matrix(object$coefficients$gamma_global[draw, ], k_domestic)) %*% W
    }
    object$coefficients$a_global <- coda::mcmc(t(temp), start = mc_start, thin = mc_thin) 
    rm(temp)
  }
  
  if (k_res > 0 | k_unres > 0) {
    object$coefficients[["c"]] <- NULL
    temp <- NULL
    if (k_unres > 0) {
      temp <- cbind(temp, object$coefficients$c_unres)
    }
    if (use_rr & k_res > 0) {
      temp <- cbind(temp, object$coefficients$c_res)
    }
    if (!is.null(temp)) {
      object$coefficients[["c"]] <- coda::mcmc(temp, start = mc_start, thin = mc_thin) 
    }
  }
  
  return(object)
}
