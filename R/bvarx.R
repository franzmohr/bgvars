.bvarx <- function(object, iterations, burnin, thin){
  #### Specifications ####
  estimation_data <- .gen_varx(object)
  y <- estimation_data$y
  x <- estimation_data$x
  tt <- dim(y)[2]
  
  k_domestic <- estimation_data$domestic$dim
  p_domestic <- estimation_data$domestic$lag
  k_foreign <- estimation_data$foreign$dim
  p_foreign <- estimation_data$foreign$lag
  
  k_x <- nrow(x)
  n_domestic <- k_domestic^2
  n_foreign <- k_domestic * k_foreign
  n_x <- k_domestic * k_x
  
  global <- !is.null(estimation_data$global)
  k_global <- 0
  n_global <- 0
  if (global) {
    k_global <- estimation_data$global$dim
    p_global <- estimation_data$global$lag
    n_global <- k_domestic * k_global
  }
  
  k_det <- estimation_data$deterministics$dim
  n_det <- k_domestic * k_det
  
  sur <- FALSE
  
  #### Priors ####
  a_mu_prior <- object$priors$coefficients$mu
  a_v_i_prior <- object$priors$coefficients$v_i
  
  sigma_df_prior <- object$priors$sigma$df
  sigma_df_post <- sigma_df_prior + tt
  sigma_scale_prior <- object$priors$sigma$scale  
  
  #### Variable selection ####
  var_select <- object$model$variable_selection$type != "none"
  a_update <- FALSE
  if (var_select) {
    var_select_type <- object$model$variable_selection$type
    a_select_draws <- iterations - burnin
    a_update <- TRUE
    a_lambda <- diag(1, n_x)
    a_prob_prior <- object$priors$variable_selection$inprior
    a_include <- object$model$variable_selection$include
    if (var_select_type == "bvs") {
      sur <- TRUE
    }
    if (var_select_type == "ssvs") {
      tau0 <- object$priors$variable_selection$tau0
      tau1 <- object$priors$variable_selection$tau1
      a_v_i_prior <- solve(diag(c(tau1^2)))
    }
    if (!is.null(object$model$variable_selection$threshold)) {
      vs_threshold <- object$model$variable_selection$threshold
      if (!is.null(object$model$variable_selection$draws)) {
        a_select_draws <- object$model$variable_selection$draws
      } else {
        a_select_draws <- floor((iterations - burnin) / 2)
      }
      if (a_select_draws > iterations - burnin) {
        stop("Number of variable selection draws may not be higher than the number of post-burn-in draws.")
      }
    }
  }
  
  #### Initial values ####
  ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
  u <- y - ols %*% x
  sigma <- tcrossprod(u) / (tt - nrow(x))
  sigma_i <- solve(sigma)
  rm(list = c("u", "ols"))
  
  if (sur) {
    z <- kronecker(t(x), diag(1, k_domestic))
    
    if (var_select) {
      z_bvs <- z
    }
  }
  
  #### Data container ####
  store <- iterations - burnin
  
  draws_a <- matrix(NA, n_x, store) 
  draws_sigma <- matrix(NA, n_domestic, store)
  draws_ll <- matrix(NA, tt, store)
  
  if (var_select) {
    draws_lambda <- matrix(NA, n_x, a_select_draws)
  }
  
  #### Gibbs sampler ####
  for (draw in 1:iterations) {
    #### Coefficients ####
    if (sur) {
      if (var_select & a_update) {
        z <- z_bvs %*% a_lambda
      }
      
      a <- post_normal_sur(y, z, sigma_i, a_mu_prior, a_v_i_prior)
      
      if (var_select) {
        if (a_update) {
          a_lambda <- bvs(y, z_bvs, a, a_lambda, sigma_i, a_prob_prior, a_include)
        }
        a <- a_lambda %*% a
      }
    } else {
      a <- post_normal(y, x, sigma_i, a_mu_prior, a_v_i_prior) 
      
      if (var_select & a_update) {
        if (var_select_type == "ssvs") {
          temp_svss <- ssvs(a, tau0, tau1, a_prob_prior, a_include)
          a_lambda <- diag(c(temp_svss$lambda))
          a_v_i_prior <- temp_svss$V_i 
        }
      }
    }
    
    #### Errors ####
    y_star <- y - matrix(a, k_domestic) %*% x
    sigma_scale_i_post <- solve(sigma_scale_prior + tcrossprod(y_star))
    sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_i_post)[,, 1], k_domestic)
    sigma <- solve(sigma_i)
    
    #### Save parameters ####
    if (draw > burnin) {
      pos_draw <- draw - burnin
      
      draws_sigma[, pos_draw] <- sigma
      draws_a[, pos_draw] <- a[1:n_x,]
      draws_ll[, pos_draw] <- loglik_normal(y_star, sigma)
      
      if (var_select) {
        if (pos_draw <= a_select_draws) {
          draws_lambda[, pos_draw] <- diag(a_lambda)
          if (pos_draw == a_select_draws & draw < iterations) {
            diag(a_lambda) <- as.numeric(rowMeans(draws_lambda) >= vs_threshold)
            if (var_select_type == "bvs") {
              z <- z_bvs %*% a_lambda
            }
            if (var_select_type == "ssvs") {
              diag(a_v_i_prior)[diag(a_lambda) == 0] <- 1 / (tau0[diag(a_lambda) == 0,]^2)
              diag(a_v_i_prior)[diag(a_lambda) == 1] <- 1 / (tau1[diag(a_lambda) == 1,]^2)
            }
            a_update <- FALSE
          }
        }
      }
      
      #### Print draws ####
      if (FALSE) {
        if (draw %in% floor(seq(from = burnin + 10, to = iterations, length.out = 100))) {
          if (pos_draw >= 10) {
            
            graphics::par(mfcol = c(2, 2))
            temp_pos <- round(.5 * pos_draw):pos_draw
            
            imp <- "s.y"
            res <- "y"
            pos_var <- k_domestic * (which(dimnames(x)[[1]] == imp) - 1) + which(dimnames(y)[[1]] == res)
            plot.ts(draws_a[pos_var, temp_pos], ylab = "contemp")
            
            res <- "y"
            if (k_det > 0) {
              pos_var <- n_x - n_det + which(dimnames(y)[[1]] == res)
              plot.ts(draws_a[pos_var, temp_pos], ylab = "deterministic")
            }
            
            imp <- "Dp"
            res <- "Dp"
            pos_var <- k_domestic * (which(dimnames(y)[[1]] == imp) - 1) + which(dimnames(y)[[1]] == res)
            plot.ts(draws_sigma[pos_var, temp_pos], ylab = "covariance") 

            
            res <- "y"
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
    pos_thin <- seq(from = a_select_draws + thin, to = store, by = thin)
    burnin <- burnin + a_select_draws
  } else {
    pos_thin <- seq(from = thin, to = store, by = thin)
  }
  
  #### Combine coefficients ####
  mc_start <- burnin + thin
  mc_thin <- thin
  
  tot_domestic <- 0
  if (p_domestic > 0) {
    tot_domestic <- n_domestic * p_domestic
    object$coefficients$a_domestic <- coda::mcmc(t(matrix(draws_a[1:tot_domestic, pos_thin], tot_domestic)), start = mc_start, thin = mc_thin) 
  }
  
  tot_foreign <- n_foreign * (p_foreign + 1)
  object$coefficients$a_foreign <- coda::mcmc(t(matrix(draws_a[tot_domestic + 1:tot_foreign, pos_thin], tot_foreign)), start = mc_start, thin = mc_thin) 
  
  tot_global <- 0
  if (global) {
    tot_global <- n_global * (p_global + 1)
    object$coefficients$a_global <- coda::mcmc(t(matrix(draws_a[tot_domestic + tot_foreign + 1:tot_global, pos_thin], tot_global)), start = mc_start, thin = mc_thin) 
  }
  
  if (k_det > 0) {
    tot_det <- k_domestic * k_det
    object$coefficients$c <- coda::mcmc(t(matrix(draws_a[tot_domestic + tot_foreign + tot_global + 1:tot_det, pos_thin], tot_det)), start = mc_start, thin = mc_thin) 
  }
  
  object$coefficients$sigma <- coda::mcmc(t(matrix(draws_sigma[, pos_thin], n_domestic)), start = mc_start, thin = mc_thin)
  
  if (var_select) {
    lambda_mean <- matrix(rowMeans(draws_lambda), k_domestic)
    dimnames(lambda_mean) <- list(dimnames(y)[[1]], dimnames(x)[[1]])
    object$coefficients$lambda_mean <- lambda_mean
  }
  
  tot_pars <- n_x
  draws_ll <- draws_ll[, pos_thin]
  LL <- sum(log(rowMeans(exp(draws_ll))))
  object$teststats <- c("loglik" = LL,
                        "AIC" = 2 * tot_pars - 2 * LL,
                        "BIC" = log(tt) * tot_pars - 2 * LL,
                        "HQ" = 2 * log(log(tt)) * tot_pars - 2 * LL)
  
  return(object)
}
