#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a BVAR model, which was produced by
#' function \code{\link{create_varxsubmodel}}.
#'
#' @param object a list of class 'varxsubmodel'.
#' @param coef a named list of prior specifications for the coefficients. It has
#' no default and must be given. It requires either \code{v_i}, the prior precision
#' of the coefficients, or \code{minnesota}; models with time varying parameters also
#' require \code{shape} and \code{rate}. Prior variances are given as precisions,
#' i.e. as their inverses. See 'Details'.
#' @param sigma a named list of prior specifications for the error term. It has no
#' default and must be given. Which elements it requires depends on the error
#' specification of the model; see 'Details'.
#' @param varsel a named list of prior specifications for the variable selection
#' algorithm. Required if the model was created with \code{varsel = "bvs"} or
#' \code{"ssvs"}, and refused if it was created with \code{varsel = "none"}.
#' See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The method has no default prior. Each argument is a named list, and
#' every value the prior needs must be given in the list it belongs to. Which
#' elements are required depends on how the sub-model was specified in
#' \code{\link{create_varxsubmodel}}: whether its parameters are time varying,
#' which error specification it has and which variable selection algorithm it
#' uses. A missing required element raises an error, and so does an element of
#' \code{coef} that the method does not recognise. The check is strict because a
#' misspelt or misplaced element, such as a shape parameter of the error term put
#' into \code{coef}, would otherwise leave the prior different from the one
#' intended without any sign of it.
#'
#' Argument \code{coef} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{a non-negative numeric specifying the prior precision of the
#'   coefficients. A value of 0 gives a flat prior. Required unless \code{minnesota} is
#'   given, and required in any case if the error covariances are estimated (error
#'   \code{"gamma+covar"} or \code{"sv+covar"}), where it is also the prior precision
#'   of the covariance coefficients.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients
#'   corresponding to deterministic terms. If it is not given, \code{v_i} is used.}
#'   \item{\code{coint_var}}{a logical specifying whether the prior mean of the first own lag of an
#'   endogenous variable should be set to 1, which is commonly used for cointegrated
#'   VAR models. If it is not given, the prior means of all lags are zero.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{coef$const = "mean"}, the mean of the respective endogenous variable is used as prior mean.
#'   If \code{coef$const = "first"}, the first values of the respective endogenous variable is used as prior mean.
#'   This can be useful for trend estimation when using a TVP model. If it is not given,
#'   the prior mean is zero.}
#'   \item{\code{minnesota}}{a named list containing the parameters for the calculation of
#'   the Minnesota prior, where the element names must be \code{kappa1}, \code{kappa2}, \code{kappa3} and \code{kappa4},
#'   of which \code{kappa1}, \code{kappa2} and \code{kappa4} are required.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}}
#'   \deqn{ \frac{\kappa_{1} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{1} \kappa_{3}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{1} \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.}
#'   \item{\code{max_var}}{accepted, but not used by this method.}
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter of the error term of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate_det}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation for coefficients, which correspond to deterministic terms. If it is not given,
#'   \code{rate} is used. Only used for models with time varying parameters.}
#'   \item{\code{omega_v}}{a positive numeric, in place of \code{shape} and \code{rate}: the
#'   variance of a normal prior on the signed standard deviation \eqn{\omega} of the state
#'   innovations, \eqn{\omega \sim N(0, V_\omega)}, which is also the prior mean of the state
#'   variance \eqn{\omega^2}. This is the non-centred parameterisation of
#'   Frühwirth-Schnatter and Wagner (2010); see Details. Only for models with time
#'   varying parameters, and not to be combined with \code{shape} and \code{rate}.}
#' }
#' If \code{minnesota} is given, the prior precisions of the coefficients are taken from it
#' and \code{v_i} and \code{v_i_det} do not enter them. The prior means are still set by
#' \code{coint_var} and \code{const}.
#'
#' Argument \code{sigma} must contain the elements that the error specification of
#' the model requires:
#' \itemize{
#'   \item error \code{"wishart"}: \code{df} and \code{scale}. A Wishart prior is refused
#'   for structural models, which need a gamma prior or stochastic volatility;
#'   \item error \code{"gamma"} or \code{"gamma+covar"}: \code{shape} and \code{rate};
#'   \item error \code{"sv"} or \code{"sv+covar"}: \code{mu}, \code{v_i}, \code{shape},
#'   \code{rate}, \code{state_variance} and \code{offset}.
#' }
#' Elements that the error specification does not use are ignored. The elements have
#' the following meaning:
#' \describe{
#'   \item{\code{df}}{a non-negative number or a character specifying the prior degrees of
#'   freedom of the inverse Wishart prior. A character is evaluated with \code{k} standing
#'   for the number of endogenous variables of the model, so that \code{"k"} or
#'   \code{"k + 3"} adapts to each sub-model, whereas a number sets the same degrees of
#'   freedom in all of them.}
#'   \item{\code{scale}}{a positive numeric specifying the prior error variance of the
#'   endogenous variables in the inverse Wishart prior.}
#'   \item{\code{shape}}{a non-negative numeric or a character specifying the prior shape
#'   parameter of the error term. A character is evaluated with \code{k} as for \code{df}.
#'   For models with stochastic volatility this prior refers to the error variance of the
#'   state equation.}
#'   \item{\code{rate}}{a positive numeric specifying the prior rate parameter of the error
#'   term. For models with stochastic volatility this prior refers to the error variance of
#'   the state equation.}
#'   \item{\code{mu}}{numeric of the prior mean of the initial state of the log-volatilities.}
#'   \item{\code{v_i}}{numeric of the prior precision of the initial state of the log-volatilities.}
#'   \item{\code{state_variance}}{numeric of the initial draw for the variance of the log-volatilities.}
#'   \item{\code{offset}}{numeric of the constant, which is added before taking the log of the squared errors.}
#' }
#'
#' For a sub-model with discounted parameters the prior is built by
#' \code{\link[bvartools]{add_priors.bvarmodel}}. It requires \code{coef$v_i},
#' \code{sigma$df} and \code{sigma$scale}, accepts only \code{v_i}, \code{v_i_det},
#' \code{v_i_alpha} and \code{const} in \code{coef} and only \code{df} and
#' \code{scale} in \code{sigma}, raises an error on any other element, and refuses
#' \code{varsel}.
#'
#' Argument \code{varsel} can contain the following elements, of which BVS requires
#' \code{inprior} and SSVS requires \code{inprior} and either \code{tau} or
#' \code{semiautomatic}:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if the variable selection algorithm
#'   should also be applied to the error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should
#'   be excluded from the variable selection algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters
#'   for the calculation of the Minnesota-like inclusion priors. See below.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second. Only used for SSVS.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables, respectively.
#'   This is the semiautomatic approach described in George et al. (2008). Only used for SSVS.}
#' }
#' 
#' If \code{varsel$minnesota} is specified, prior
#' inclusion probabilities are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_{1}}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_{2}}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_{3}}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3},
#' \eqn{\kappa_4} as the first, second, third and forth element in
#' \code{varsel$minnesota}, respectively.
#' 
#' The state variances of a model with time varying parameters take one of two
#' parameterisations, and \code{coef} names which by naming its parameters. With
#' \code{shape} and \code{rate} the precisions of the state equation have a gamma
#' prior, so the variances have an inverse gamma one, which gives zero density to
#' a variance of zero: the model cannot conclude that a coefficient is constant,
#' because the prior has already ruled it out. With \code{omega_v} the random walk
#' \eqn{x_t = x_{t-1} + v_t}, \eqn{v_t \sim N(0, \omega^2)}, is written
#' \eqn{x_t = x_0 + \omega \tilde{x}_t} with \eqn{\tilde{x}_t} a standard random
#' walk, so \eqn{\omega} enters the measurement equation as a regression
#' coefficient and a normal prior \eqn{\omega \sim N(0, V_\omega)} is conjugate.
#' A coefficient that does not move, \eqn{\omega = 0}, is then an interior point
#' of the prior rather than the boundary of a variance, and the posterior carries
#' the Savage-Dickey density ratio that tests it: the draws of \eqn{\omega} come
#' back in \code{posterior$a$omega} and the log density at zero of its conditional
#' posterior in \code{posterior$a$omega_log_zero}, one column per state, with
#' \code{omega_log_zero_joint} for the block as a whole. This is the non-centred
#' parameterisation of Frühwirth-Schnatter and Wagner (2010), and the test is that
#' of Chan (2018).
#'
#' \eqn{V_\omega} is also the prior mean of the state variance, the implied prior
#' on \eqn{\omega^2} being \eqn{Gamma(1/2, 1 / (2 V_\omega))}, so the two
#' parameterisations are given the same prior mean by setting \code{omega_v} to
#' the \code{rate / (shape - 1)} of the other. What differs is the shape: the
#' non-centred prior puts far more mass near zero. They cannot be combined, and
#' the covariance block follows whichever the coefficients take.
#'
#' @return An object of class 'varxsubmodel'.
#' 
#' @references
#' 
#' Chan, J. C. C. (2018). Specification tests for time-varying parameter models with stochastic
#' volatility. \emph{Econometric Reviews, 37}(8), 807--823. \doi{10.1080/07474938.2016.1167948}
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#'
#' Frühwirth-Schnatter, S., & Wagner, H. (2010). Stochastic model specification search for
#' Gaussian and partial non-Gaussian state space models. \emph{Journal of Econometrics, 154}(1),
#' 85--100. \doi{10.1016/j.jeconom.2009.07.003}
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples
#' 
#' # Load data
#' data("gvar2019")
#' global_data <- gvar2019[["global_data"]]
#' submodel_data <- gvar2019[["submodel_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#' 
#' # Create global model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' 
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2013:2016)
#' 
#' # Create sub-model
#' object <- create_varxsubmodel(object,
#'                               submodel = "AT",
#'                               endogen = c("y","Dp", "r"), p_endogen = 1,
#'                               exogen = c("y", "Dp"), p_exogen = 2,
#'                               global = "poil", s = 1,
#'                               deterministic = "const", seasonal = FALSE,
#'                               structural = FALSE, tvp = FALSE,
#'                               error = "wishart", varsel = "none",
#'                               iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # The previous function returns a model list. Extract the first model to proceed.
#' object <- object[[1]]
#' 
#' # Add priors
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 1e-8))
#'                               
#' 
#' 
#' @export
#' @method add_priors varxsubmodel
add_priors.varxsubmodel <- function(object,
                                    coef,
                                    sigma,
                                    varsel = NULL,
                                    ...){
  
  # The discounted model takes a matrix normal prior over the coefficient matrix
  # rather than a normal one over its vectorisation. A sub-model's prior is the
  # same object as a bvarmodel's, so the method there builds it.
  if (bvartools::is_discount_model(object)) {
    return(NextMethod())
  }

  # Input checks
  ## Coefficient priors ----
  if (!is.null(coef)) {
    
    .add_priors_check_coef(object, coef)
    
    if (!is.null(coef[["const"]])) {
      if ("character" %in% class(coef[["const"]])) {
        if (!coef[["const"]] %in% c("first", "mean")) {
          stop("Invalid specificatin of coef$const.")
        }
      }
    }
    
    if (is.null(coef[["v_i_det"]])) {
      coef[["v_i_det"]] <- coef[["v_i"]]
    }
    
  }
  
  ## Sigma priors ----
  error_prior <- .add_priors_check_sigma(object, sigma)
  
  # Check Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  # Check coint VAR ----
  coint_var <- FALSE # Cointegrated VAR?
  if (!is.null(coef[["coint_var"]])) {
    if (coef[["coint_var"]]) {
      coint_var <- TRUE 
    }
  }
  
  if (!is.null(varsel) & object[["model"]][["varsel"]] == "none") {
    stop("Argument 'varsel' was provided, but according to argument 'object$model$varsel' no variable selection algorithm should be used.")
  }
  
  # Check SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (object[["model"]][["varsel"]] == "ssvs") {
    
    .add_priors_check_ssvs(object, varsel)
    
    if (!is.null(varsel[["covar"]])) {
      if (varsel[["covar"]]) {
        use_ssvs_error <- TRUE 
      }
    }
    
    if (is.null(varsel[["exclude_det"]])) {
      varsel[["exclude_det"]] <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approach of 
    # George et al. (2008) should be used
    if (!is.null(varsel[["semiautomatic"]])) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
  }
  
  # BVS ----
  use_bvs <- FALSE
  use_bvs_error <- FALSE
  if (object[["model"]][["varsel"]] == "bvs") {
    use_bvs <- TRUE
    
    .add_priors_check_bvs(object, varsel)
    
    if (is.null(varsel[["exclude_det"]])) {
      varsel[["exclude_det"]] <- FALSE
    }
    if (!is.null(varsel[["covar"]])) {
      if (varsel[["covar"]]) {
        use_bvs_error <- TRUE 
      }
    }
    .add_priors_warn_bvs_flat(coef, minnesota, varsel[["exclude_det"]])
  }
  
  if (use_ssvs & use_bvs) {
    stop("SSVS and BVS cannot be applied at the same time.")
  }
  
  if (error_prior == "wishart" & (use_ssvs_error | use_bvs_error)) {
    stop("Wishart prior not allowed when BVS or SSVS are applied to covariance matrix.")
  }
  
  varsel_covar <- use_ssvs_error | use_bvs_error
  
  # Generate priors ----
  
  # Get model specs to obtain total number of coeffs
  k <- object[["model"]][["k_endogen"]]
  p <- object[["model"]][["p_endogen"]]
  
  if (k == 1 & (use_ssvs_error | use_bvs_error)) {
    stop("BVS or SSVS cannot be applied to covariance matrix when there is only one endogenous variable.")
  } 
  
  # 'p_exogen' and 's_global' count blocks of regressors from the contemporaneous
  # term onwards, so they already are the number of blocks and are used as they
  # are. A value of zero means that the variables contribute no regressors.
  if (object[["model"]][["k_exogen"]] > 0) {
    k_exogen <- object[["model"]][["k_exogen"]]
    p_exogen <- object[["model"]][["p_exogen"]]
  } else {
    k_exogen <- 0
    p_exogen <- 0
  }
  
  if (object[["model"]][["m_global"]] > 0) {
    m <- object[["model"]][["m_global"]]
    s <- object[["model"]][["s_global"]]
  } else {
    s <- 0
    m <- 0
  }
  
  # Total # of non-deterministic coefficients
  n_endogen <- k * (k * p)
  n_exogen <- k * (k_exogen * p_exogen)
  n_global <- k * (m * s)
  
  # Add number of non-cointegration deterministic terms
  n_det <- 0
  if (object[["model"]][["n"]] > 0){
    n_det <- object[["model"]][["n"]] * k
  }
  
  tot_par <- n_endogen + n_exogen + n_global + n_det
  
  covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  structural <- object[["model"]][["structural"]]
  if (covar & structural) {
    stop("Error covariances and structural coefficients cannot be estimated at the same time.")
  }
  
  # A constant coefficient sampler selects over one set of coefficients or over
  # both: it reads a single selection scheme for the whole model, so a
  # covariance block it is given goes into the selection with the rest. Only the
  # time varying samplers take the covariance block's scheme separately, which
  # is why the same call is allowed there. Left to run, this combination fails
  # inside the sampler on a prior it was never given.
  if ((use_ssvs | use_bvs) & covar & !varsel_covar &
      !object[["model"]][["tvp"]] & k > 1) {
    stop("Variable selection cannot be restricted to the coefficients when the ",
         "model has an error covariance block and constant coefficients: this ",
         "sampler applies one selection scheme to both. Set 'varsel$covar' to ",
         "TRUE to select over the covariances as well, drop the covariances with ",
         "an 'error' of \"gamma\" or \"sv\", or use a time varying model, where ",
         "the two blocks can differ.")
  }
  
  n_struct <- 0
  if (structural & k > 1) {
    n_struct <- (k - 1) * k / 2
    tot_par <- tot_par + n_struct
  }
  
  # Additional input check
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    if (tot_par != ncol(object[["data"]][["train"]][["z"]])) {
      stop("Model specifications are not consistent with data matrix 'object$data$train$z'.")
    } 
  }
  
  sv <- error_prior == "sv"
  
  # Priors ----
  ## Coefficients ----
  if (tot_par > 0) {
    
    ### Prior means ----
    
    mu <- matrix(rep(0, tot_par - n_struct), k)
    
    # Add 1 to first own lags for cointegrated VAR model
    if (coint_var & p > 0) {
      mu[1:k, 1:k] <- diag(1, k)
    }
    
    # Prior for intercept terms
    if (n_det > 0) {
      
      if (!is.null(coef[["const"]]))  {
        
        pos <- which(dimnames(object[["data"]][["train"]][["x"]])[[2]] == "const")
        
        if (length(pos) == 1) {
          if ("character" %in% class(coef[["const"]])) {
            if (coef[["const"]] == "first") {
              mu[, pos] <- object[["data"]][["train"]][["y"]][1, ]
            }
            if (coef[["const"]] == "mean") {
              mu[, pos] <- colMeans(object[["data"]][["train"]][["y"]])
            }
          }
          if ("numeric" %in% class(coef[["const"]])) {
            if (length(coef[["const"]]) == 1 | length(coef[["const"]]) == k) {
              mu[, pos] <- coef[["const"]]  
            } else {
              stop("When a numeric is provided in argument 'coef$const', it must be either a single number or a vector of the same length as the number of endogenous varibles in the model.")
            }
          } 
        }
      }
    }
    
    mu <- matrix(mu)
    
    if (structural) {
      mu <- rbind(mu, matrix(0, n_struct))
    }
    
    object[["priors"]][["a"]] <- list(type = "normal",
                                      mu = mu)
    
    ### Prior covariances ----
    
    if (minnesota) {
      #### Minnesota prior ----
      minn <- minnesota_prior(object = object,
                              kappa1 = coef[["minnesota"]][["kappa1"]],
                              kappa2 = coef[["minnesota"]][["kappa2"]],
                              kappa3 = coef[["minnesota"]][["kappa3"]],
                              kappa4 = coef[["minnesota"]][["kappa4"]],
                              max_var = NULL,
                              coint_var = FALSE,
                              sigma = "AR")
      
      object[["priors"]][["a"]][["v_inv"]] <- minn[["v_inv"]]
    }
    
    #### SSVS prior ----
    if (use_ssvs) {
      
      if (object[["model"]][["tvp"]]) {
        stop("SSVS is not supported for TVP models.")
      }
      
      if (sv) {
        stop("Not allowed to use SSVS with stochastic volatility models.")
      }
      
      ssvs_temp <- ssvs_prior(object, tau = varsel[["tau"]], semiautomatic = varsel[["semiautomatic"]])
      temp <- inclusion_prior(object, prob = varsel[["inprior"]], exclude_deterministics = varsel[["exclude_det"]],
                              minnesota_like = !is.null(varsel[["minnesota"]]),
                              kappa1 = varsel[["minnesota"]][1],
                              kappa2 = varsel[["minnesota"]][2],
                              kappa3 = varsel[["minnesota"]][3],
                              kappa4 = varsel[["minnesota"]][4])
      object[["priors"]][["a"]][["v_inv"]] <- diag(1 / ssvs_temp[["tau1"]][, 1]^2, tot_par)
      object[["priors"]][["a"]][["inprior"]] <- temp[["prior"]]
      object[["priors"]][["a"]][["include"]] <- temp[["include"]]
      object[["priors"]][["a"]][["tau0"]] <- ssvs_temp[["tau0"]]
      object[["priors"]][["a"]][["tau1"]] <- ssvs_temp[["tau1"]]
      rm(temp)
      rm(ssvs_temp)
    }
    
    #### Regular prior ----
    if (!minnesota & !use_ssvs) {
      v_i <- diag(coef[["v_i"]], tot_par)
      # Add priors for deterministic terms if they were specified
      if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
        diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
      }
      object[["priors"]][["a"]][["v_inv"]] <- v_i
    }
    
    #### BVS prior ----
    if (use_bvs) {
      temp <- inclusion_prior(object, prob = varsel[["inprior"]],
                              exclude_deterministics = varsel[["exclude_det"]],
                              minnesota_like = !is.null(varsel[["minnesota"]]),
                              kappa1 = varsel[["minnesota"]][1],
                              kappa2 = varsel[["minnesota"]][2],
                              kappa3 = varsel[["minnesota"]][3],
                              kappa4 = varsel[["minnesota"]][4])
      object[["priors"]][["a"]][["inprior"]] <- temp[["prior"]]
      object[["priors"]][["a"]][["include"]] <- temp[["include"]]
      rm(temp)
    }
    
    ### TVP prior ----
    if (object[["model"]][["tvp"]]) {
      if (!is.null(coef[["omega_v"]])) {
        # The non-centred parameterisation; see Details, and the VEC method,
        # whose loadings take a scale of their own that a VARX has no
        # counterpart for.
        object[["priors"]][["a"]][["omega_v"]] <- matrix(coef[["omega_v"]], tot_par)
      } else {
        object[["priors"]][["a"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
        object[["priors"]][["a"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
        if (n_det > 0 & !is.null(coef[["rate_det"]])) {
          object[["priors"]][["a"]][["rate"]][tot_par - n_struct - n_det + 1:n_det, ] <- coef[["rate_det"]]
        }
      }
    }
  }
  
  ## Covar priors ----
  
  if (!structural & covar & k > 1) {
    
    n_covar <- k * (k - 1) / 2
    object[["priors"]][["psi"]][["type"]] <- "normal"
    object[["priors"]][["psi"]][["mu"]] <- matrix(0, n_covar)
    object[["priors"]][["psi"]][["v_inv"]] <- diag(coef[["v_i"]], n_covar)
    if (object[["model"]][["tvp"]]) {
      if (!is.null(coef[["omega_v"]])) {
        object[["priors"]][["psi"]][["omega_v"]] <- matrix(coef[["omega_v"]], n_covar)
      } else {
        object[["priors"]][["psi"]][["shape"]] <- matrix(coef[["shape"]], n_covar)
        object[["priors"]][["psi"]][["rate"]] <- matrix(coef[["rate"]], n_covar)
      }
    }
    
    # Variable selection
    object[["priors"]][["psi"]][["varsel"]] <- "none"
    
    # SSVS priors
    if (use_ssvs_error) {
      object[["priors"]][["psi"]][["varsel"]] <- "ssvs"
      object[["priors"]][["psi"]][["inprior"]] <- matrix(varsel[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["include"]] <- matrix(1:n_covar)
      object[["priors"]][["psi"]][["tau0"]] <- matrix(varsel[["tau"]][1], n_covar)
      object[["priors"]][["psi"]][["tau1"]] <- matrix(varsel[["tau"]][2], n_covar)
    }
    
    # BVS priors
    if (use_bvs_error) {
      object[["priors"]][["psi"]][["varsel"]] <- "bvs"
      object[["priors"]][["psi"]][["inprior"]] <- matrix(varsel[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["include"]] <- matrix(1:n_covar)
    }
  }
  
  ## Error term ----
  if (error_prior == "sv") {
    
    object <- .add_priors_sv_helper(object, sigma, k)
    
  } else {
    
    if (error_prior == "wishart") {
      object[["priors"]][["u_sigma"]][["type"]] <- "wishart"
      help_df <- sigma[["df"]]
      object[["priors"]][["u_sigma"]][["df"]] <- NA_real_
      object[["priors"]][["u_sigma"]][["scale"]] = diag(sigma[["scale"]], k)
    }
    
    if (error_prior == "gamma") {
      object[["priors"]][["u_sigma"]][["type"]] <- "gamma"
      help_df <- sigma[["shape"]]
      object[["priors"]][["u_sigma"]][["shape"]] <- NA_real_
      object[["priors"]][["u_sigma"]][["rate"]] = matrix(sigma[["rate"]], k)
    }
    
    if (minnesota & !is.null(object[["data"]][["train"]][["x"]])) {
      # Store LS estimate of variance covariance matrix for analytical solution
      object[["priors"]][["u_sigma"]][["u_sigma_inv"]] = minn[["sigma_inv"]]
    }
    
    if ("character" %in% class(help_df)) {
      if (grepl("k", help_df)) {
        # Transform character specification to expression and evaluate
        help_df <- eval(parse(text = help_df))
      } else {
        stop("Use no other letter than 'k' in 'sigma$df' to indicate the number of endogenous variables.")
      }
    }
    
    if (help_df < 0) {
      stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
    }
    
    if (error_prior == "wishart") {
      object[["priors"]][["u_sigma"]][["df"]] <- as.integer(help_df)
    }
    if (error_prior == "gamma") {
      object[["priors"]][["u_sigma"]][["shape"]] <- matrix(help_df, k)
    } 
  }
  
  return(object)
}