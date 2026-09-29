#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{create_vecxsubmodel}}.
#'
#' @param object a list of class 'vecxsubmodel'.
#' @param coef a named list of prior specifications for the coefficients that do not
#' determine the cointegration space, i.e. the loadings, the coefficients of the
#' differenced regressors and those of the unrestricted deterministic terms. It has no
#' default and must be given. It requires either \code{v_i}, the prior precision of the
#' coefficients, or \code{minnesota}; models with time varying parameters also require
#' \code{shape} and \code{rate}. Prior variances are given as precisions, i.e. as their
#' inverses. See 'Details'.
#' @param coint a named list of prior specifications for the cointegration space. It has
#' no default and must be given, except for sub-models with discounted parameters, which
#' take none. See 'Details'.
#' @param sigma a named list of prior specifications for the error term. It has no
#' default and must be given. Which elements it requires depends on the error
#' specification of the model; see 'Details'.
#' @param varsel a named list of prior specifications for the variable selection
#' algorithm. Required if the model was created with \code{varsel = "bvs"}, and refused
#' if it was created with \code{varsel = "none"}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The method has no default prior. Each argument is a named list, and
#' every value the prior needs must be given in the list it belongs to. Which
#' elements are required depends on how the sub-model was specified in
#' \code{\link{create_vecxsubmodel}}: whether its parameters are time varying,
#' which error specification it has and whether it uses variable selection. A
#' missing required element raises an error, and so does an element of \code{coef}
#' or \code{coint} that the method does not recognise. The check is strict because a
#' misspelt or misplaced element, such as a precision of the cointegration space put
#' into \code{coef}, would otherwise leave the prior different from the one intended
#' without any sign of it.
#'
#' Argument \code{coef} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{a non-negative numeric specifying the prior precision of the
#'   coefficients. A value of 0 gives a flat prior. Required unless \code{minnesota} is
#'   given, and required in any case if the error covariances are estimated (error
#'   \code{"gamma+covar"} or \code{"sv+covar"}), where it is also the prior precision
#'   of the covariance coefficients. For models with time varying parameters and a
#'   positive rank the prior precision of the loadings is \code{v_i} divided by
#'   \eqn{1 - \rho^2}, where \eqn{\rho} is \code{coint$rho}.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients
#'   corresponding to unrestricted deterministic terms. If it is not given, \code{v_i}
#'   is used.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{coef$const = "mean"}, the mean of the respective endogenous variable is used as prior mean.
#'   If \code{coef$const = "first"}, the first values of the respective endogenous variable is used as prior mean.
#'   If it is not given, the prior mean is zero.}
#'   \item{\code{minnesota}}{a named list containing the parameters for the calculation of
#'   the Minnesota prior, where the element names must be \code{kappa1}, \code{kappa2}, \code{kappa3} and \code{kappa4},
#'   of which \code{kappa1}, \code{kappa2} and \code{kappa4} are required.
#'   The prior precisions of the coefficients are then taken from it, and \code{v_i} and
#'   \code{v_i_det} do not enter them.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{1} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{1} \kappa_{3}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{1} \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   The function only provides priors for the non-cointegration part of the model. However,
#'   the residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}, \code{coint_var}}{accepted, because the VARX method takes them,
#'   but not used by this method.}
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter of the error term of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate_det}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation for coefficients, which correspond to deterministic terms. If it is not given,
#'   \code{rate} is used. Only used for models with time varying parameters.}
#'   \item{\code{rate_alpha}}{a numeric specifying the prior rate parameter of the error term of
#'   the state equation for the loadings. If it is not given, \code{rate} is used. Only used for
#'   models with time varying parameters and a positive rank. The loadings multiply levels, so a
#'   drift that is negligible for a coefficient on a differenced regressor moves the fitted
#'   value of an equation by much more; a rate several orders of magnitude below \code{rate}
#'   keeps their drift from absorbing the residuals while the other coefficients vary. It does
#'   not reach the steps of the cointegration vectors, which can absorb the residuals as well;
#'   see section 'Prior on the cointegration space'.}
#'   \item{\code{omega_v}}{a positive numeric, in place of \code{shape} and \code{rate}: the
#'   variance of a normal prior on the signed standard deviation \eqn{\omega} of the state
#'   innovations, \eqn{\omega \sim N(0, V_\omega)}, which is also the prior mean of the state
#'   variance \eqn{\omega^2}. This is the non-centred parameterisation of
#'   Frühwirth-Schnatter and Wagner (2010); see Details. Only for models with time
#'   varying parameters, and not to be combined with \code{shape} and \code{rate}.}
#'   \item{\code{omega_v_alpha}}{a positive numeric, the \code{omega_v} of the loadings, as
#'   \code{rate_alpha} is the \code{rate} of the loadings. If it is not given, \code{omega_v} is
#'   used. Only used with \code{omega_v} and a positive rank.}
#' }
#' 
#' Argument \code{coint} specifies the prior on the cointegration space. It is
#' checked and built by \code{\link[bvartools]{cointspace_prior}}, the function
#' \pkg{bvartools} uses for its own VEC models, and its section 'Prior on the
#' cointegration space', shown below, says which elements are required. The error
#' correction term \code{data$train$w} of a sub-model holds its endogenous, weakly
#' exogenous and global variables and its restricted deterministic terms, and the
#' prior has one row and column for each of them. Unrecognised elements of
#' \code{coint} raise an error.
#'
#' Variable selection is available as \code{"bvs"}. The loadings are never subject to it,
#' so the error correction term is always in the model, and one inclusion parameter is drawn
#' per coefficient rather than per coefficient and period. SSVS is not implemented for error
#' correction sub-models and is refused when the priors are added.
#' 
#' Argument \code{sigma} must contain the elements that the error specification of
#' the model requires:
#' \itemize{
#'   \item error \code{"wishart"}: \code{df} and \code{scale}. A Wishart prior is refused
#'   for structural models, which need a gamma prior or stochastic volatility;
#'   \item error \code{"gamma"} or \code{"gamma+covar"}: \code{shape} and \code{rate};
#'   \item error \code{"sv"} or \code{"sv+covar"}: \code{mu}, \code{v_i}, \code{shape},
#'   \code{rate}, \code{state_variance} and \code{offset}; with \code{tvp = TRUE},
#'   \code{omega_v} may take the place of \code{shape} and \code{rate}.
#' }
#' Elements that the error specification does not use are ignored. The elements have
#' the following meaning:
#' \describe{
#'   \item{\code{df}}{a non-negative number or a character specifying the prior degrees of
#'   freedom of the inverse Wishart prior. A character is evaluated with \code{k} standing
#'   for the number of endogenous variables of the model, so that \code{"k"} or
#'   \code{"k + 3"} adapts to each sub-model, whereas a number sets the same degrees of
#'   freedom in all of them. The value is stored as given. The sampler adds the rank
#'   \eqn{r} of the cointegration matrix to the posterior degrees of freedom, as the prior
#'   of the loadings requires.}
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
#'   \item{\code{omega_v}}{a positive numeric, in place of \code{shape} and \code{rate}: the
#'   variance of a normal prior on the signed standard deviation of the log-volatility
#'   innovations, as \code{coef$omega_v} is for the coefficients. Only for sub-models with
#'   time varying parameters and stochastic volatility, and never beside \code{shape} and
#'   \code{rate}. Under it the log-volatility does not move exactly when the standard
#'   deviation is zero, which is a point inside the prior, so
#'   \code{\link[bvartools]{time_variation_test}} can report a Bayes factor for that block
#'   as well as for the coefficients.}
#' }
#' 
#' For a sub-model with discounted parameters the prior is built by
#' \code{\link[bvartools]{add_priors.bvecmodel}}. Such a model conditions on a fixed
#' cointegration space rather than drawing one, so it takes no \code{coint}, and
#' giving one raises an error. It requires \code{coef$v_i}, \code{sigma$df} and
#' \code{sigma$scale}, accepts only \code{v_i}, \code{v_i_det}, \code{v_i_alpha} and
#' \code{const} in \code{coef} and only \code{df} and \code{scale} in \code{sigma},
#' raises an error on any other element, and refuses \code{varsel}.
#' 
#' Argument \code{varsel} can contain the following elements, of which BVS requires
#' \code{inprior}:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if the variable selection algorithm
#'   should also be applied to the error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should
#'   be excluded from the variable selection algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters
#'   for the calculation of the Minnesota-like inclusion priors. See below.}
#' }
#' 
#' If \code{varsel$minnesota} is specified, prior
#' inclusion probabilities are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for exogenous variables, \cr
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
#' @inheritSection bvartools::cointspace_prior Prior on the cointegration space
#'
#' @return An object of class 'vecxsubmodel'.
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
#' Johansen, S. (1995). \emph{Likelihood-based inference in cointegrated vector
#' autoregressive models}. Oxford: Oxford University Press.
#'
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples
#' 
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("EA", "US"))
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
#' object <- add_submodels(object,
#'                         p_endogen = 1, p_exogen = 2,
#'                         global = "poil",
#'                         s = 2,
#'                         r = 1,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' object <- add_priors(object[["submodels"]][["EA"]][[1]],
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' 
#' @export
#' @method add_priors vecxsubmodel
add_priors.vecxsubmodel <- function(object,
                                    coef,
                                    coint,
                                    sigma,
                                    varsel = NULL,
                                    ...){
  
  # The discounted model takes a matrix normal prior over the coefficient matrix
  # rather than a normal one over its vectorisation, and no cointegration space
  # prior at all: it conditions on a fixed space rather than drawing one. A
  # sub-model's prior is the same object as a bvecmodel's, so the method there
  # builds it, and this one gets out of the way before `coint`, which is never
  # given for it, is forced.
  if (bvartools::is_discount_model(object)) {
    return(NextMethod())
  }

  # Input checks ----
  
  ## coefficients ----
  if (!is.null(coef)) {
    .add_priors_check_coef(object, coef)
  }
  
  ## cointegration ----
  # Checked and built by bvartools, whose VEC models get the same prior. Its
  # dimensions come from the error correction term data$train$w, which for a
  # sub-model also holds the weakly exogenous and global variables. It is stored
  # further down, with the others.
  beta_prior <- bvartools::cointspace_prior(object, coint)
  
  ## sigma ----
  error_prior <- .add_priors_check_sigma(object, sigma)
  
  ## Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  if (!is.null(varsel) & object[["model"]][["varsel"]] == "none") {
    stop("Argument 'varsel' was provided, but according to argument 'object$model$varsel' no variable selection algorithm should be used.")
  }
  
  ## SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (object[["model"]][["varsel"]] == "ssvs") {
    stop("Stochastic search variable selection is not implemented for error ",
         "correction sub-models. Consider using 'bvs' instead.")
    #.add_priors_check_ssvs(object, varsel)
    
    if (!is.null(varsel[["covar"]])) {
      if (varsel[["covar"]]) {
        use_ssvs_error <- TRUE 
      }
    }
    
    if (is.null(varsel[["exclude_det"]])) {
      varsel[["exclude_det"]] <- FALSE
    }
    # In case varsel is specified, check if the semi-automatic approach of 
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
  
  ## BVS ----
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
  
  # A constant coefficient sampler selects over one set of coefficients or over
  # both: it reads a single selection scheme for the whole model, so a
  # covariance block it is given goes into the selection with the rest. Only the
  # time varying samplers take the covariance block's scheme separately, which
  # is why the same call is allowed there. Left to run, this combination fails
  # inside the sampler on a prior it was never given.
  if ((use_ssvs | use_bvs) &
      object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") &
      !varsel_covar & !object[["model"]][["tvp"]] &
      object[["model"]][["k_endogen"]] > 1) {
    stop("Variable selection cannot be restricted to the coefficients when the ",
         "model has an error covariance block and constant coefficients: this ",
         "sampler applies one selection scheme to both. Set 'varsel$covar' to ",
         "TRUE to select over the covariances as well, drop the covariances with ",
         "an 'error' of \"gamma\" or \"sv\", or use a time varying model, where ",
         "the two blocks can differ.")
  }
  
  # Generate priors ----
  
  # Get model specs to obtain total number of coeffs
  k_endogen <- object[["model"]][["k_endogen"]]
  p_endogen <- object[["model"]][["p_endogen"]]
  k_exogen <- object[["model"]][["k_exogen"]]
  p_exogen <- object[["model"]][["p_exogen"]]
  m_global <- object[["model"]][["m_global"]]
  s_global <- object[["model"]][["s_global"]]
  r <- object[["model"]][["rank"]]
  
  if (k_endogen == 1 & (use_ssvs_error | use_bvs_error)) {
    stop("BVS or SSVS cannot be applied to covarianc matrix when there is only one endogenous variable.")
  } 
  
  use_global <- m_global > 0
  
  # Substract lag from domestic model for VEC
  p_endogen <- p_endogen - 1
  
  # The weakly exogenous and the global variables carry the lag orders of the
  # model in levels as well, counted from the contemporaneous term onwards, so
  # the model holds one block of their differences less.
  p_exogen <- max(p_exogen - 1, 0)
  s_global <- max(s_global - 1, 0)
  
  # Total # of non-deterministic coefficients
  n_alpha <- k_endogen * r
  n_gamma_endogen <- k_endogen * (k_endogen * p_endogen)
  n_gamma_exogen <- k_endogen * (k_exogen * p_exogen)
  n_upsilon <- k_endogen * (m_global * s_global)
  
  # Add number of unrestricted deterministic terms
  n_det <- k_endogen * object[["model"]][["n"]]
  
  tot_par <- n_alpha + n_gamma_endogen + n_gamma_exogen + n_upsilon + n_det
  
  covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  structural <- object[["model"]][["structural"]]
  if (covar & structural) {
    stop("Error covariances and structural coefficients cannot be estimated at the same time.")
  }
  sv <- object[["model"]][["error"]] %in% c("sv", "sv+covar")
  n_struct <- 0
  n_z <- NCOL(object[["data"]][["train"]][["z"]])
  if (object[["model"]][["rank"]] > 0) {
    n_w <- NCOL(object[["data"]][["train"]][["w"]])
  }
  
  if (structural & k_endogen > 1) {
    n_struct <- (k_endogen - 1) * k_endogen / 2
    tot_par <- tot_par + n_struct
  }
  
  # Additional input check
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    if (tot_par != ncol(object[["data"]][["train"]][["z"]])) {
      stop("Model specifications are not consistent with data matrix 'object$data$train$z'.")
    } 
  }
  
  #### Cointegration ----
  if (r > 0) {
    object[["priors"]][["beta"]] <- beta_prior
  }
  
  # Non-cointegration ----
  # Generate prior matrices ----
  if (tot_par > 0) {
    
    # Zero prior means
    mu <- matrix(rep(0, tot_par - n_struct), k_endogen)
    
    # Prior for intercept terms
    if (n_det > 0) {
      
      if (!is.null(coef[["const"]]))  {
        
        # The column of the constant among the regressors of the
        # non-cointegration term. The columns of 'mu' start with the r of the
        # loadings, so the constant is column r + pos. 'data$x' does not exist,
        # and looking the constant up there left every specification of
        # coef$const without effect.
        pos <- which(dimnames(object[["data"]][["train"]][["x"]])[[2]] == "const")

        if (length(pos) == 1) {
          if ("character" %in% class(coef[["const"]])) {
            if (coef[["const"]] == "first") {
              mu[, r + pos] <- object[["data"]][["train"]][["y"]][1,]
            }
            if (coef[["const"]] == "mean") {
              mu[, r + pos] <- colMeans(object[["data"]][["train"]][["y"]])
            }
          }
          if ("numeric" %in% class(coef[["const"]])) {
            if (length(coef[["const"]]) == 1 | length(coef[["const"]]) == k_endogen) {
              mu[, r + pos] <- coef[["const"]]
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
    
    # Prior covariances
    if (minnesota) {
      # Minnesota prior ----
      minn <- minnesota_prior(object = object,
                              kappa1 = coef[["minnesota"]][["kappa1"]],
                              kappa2 = coef[["minnesota"]][["kappa2"]],
                              kappa3 = coef[["minnesota"]][["kappa3"]],
                              kappa4 = coef[["minnesota"]][["kappa4"]],
                              max_var = NULL,
                              coint_var = FALSE,
                              sigma = "AR")
      
      # minnesota_prior.vecxsubmodel() calls the precision 'v_i', as the method
      # of bvartools for VEC models does. Reading 'v_inv' stored NULL, so a
      # Minnesota prior left the sampler without a prior precision.
      object[["priors"]][["a"]][["v_inv"]] <- minn[["v_i"]]
    }
    
    # SSVS prior ----
    # Unreachable while the specification above refuses SSVS, and left in place
    # so that the branch is ready when it is implemented.
    if (use_ssvs) {
      stop("Stochastic search variable selection is not implemented for error ",
           "correction sub-models. Consider using 'bvs' instead.")
      if (object[["model"]][["tvp"]]) {
        stop("SSVS is not supported for TVP models.")
      }
      
      if (sv) {
        stop("Not allowed to use SSVS with stochastic volatility models.")
      }
      
      ssvs_temp <- ssvs_prior(object, tau = varsel[["tau"]], semiautomatic = varsel[["semiautomatic"]])
      temp <- inclusion_prior(object, prob = varsel[["inprior"]], exclude_deterministics = varsel[["exclude_det"]],
                              minnesota_like = !is.null(varsel[["minnesota"]]), kappa = varsel[["minnesota"]])
      object[["priors"]][["a"]][["v_inv"]] <- diag(1 / ssvs_temp[["tau1"]][, 1]^2, tot_par)
      object[["priors"]][["a"]][["inprior"]] <- temp[["prior"]]
      object[["priors"]][["a"]][["include"]] <- temp[["include"]]
      object[["priors"]][["a"]][["tau0"]] <- ssvs_temp[["tau0"]]
      object[["priors"]][["a"]][["tau1"]] <- ssvs_temp[["tau1"]]
      rm(temp)
      rm(ssvs_temp)
    }
    
    # Regular prior ----
    if (!minnesota & !use_ssvs) {
      
      if (object[["model"]][["tvp"]]) {

        # For a TVP model this is the prior precision of the state before the
        # sample rather than of a constant coefficient.
        v_i <- diag(coef[["v_i"]], tot_par)

        # The loadings carry the compensating scale of the cointegration space.
        # beta_t has stationary variance 1 / (1 - rho^2), which for a rho just
        # below one is large, and only the product alpha beta' is identified --
        # so alpha's prior variance is shrunk by the same factor, leaving the
        # product on the scale coef$v_i asks for. Koop et al. (2011) scale the
        # loadings' variance by 1 - rho^2 relative to the other coefficients,
        # which is why coef$v_i enters here: a precision of 1 / (1 - rho^2)
        # alone put the product at unit variance whatever coef$v_i was.
        if (r > 0) {
          diag(v_i)[1:n_alpha] <- coef[["v_i"]] / (1 - coint[["rho"]] * coint[["rho"]])
        }
        if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
        }
        if (!is.null(coef[["omega_v"]])) {
          # The non-centred parameterisation: a normal prior on the signed
          # standard deviation of the state innovations instead of the gamma one
          # on their precision. Written this way a coefficient that does not move
          # is an interior point of the prior rather than the boundary of a
          # variance, which is what carries the test of time variation; see
          # Details. The two are mutually exclusive, and .add_priors_check_coef()
          # has already refused a 'coef' that gives both.
          object[["priors"]][["a"]][["omega_v"]] <- matrix(coef[["omega_v"]], tot_par)
          # The loadings come first among the coefficients, and multiply levels
          # rather than differences, so they take a prior of their own here as
          # they do through 'rate_alpha' in the centred parameterisation.
          if (r > 0 & !is.null(coef[["omega_v_alpha"]])) {
            object[["priors"]][["a"]][["omega_v"]][1:n_alpha, ] <- coef[["omega_v_alpha"]]
          }
        } else {
          object[["priors"]][["a"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
          object[["priors"]][["a"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
          if (n_det > 0 & !is.null(coef[["rate_det"]])) {
            object[["priors"]][["a"]][["rate"]][tot_par - n_struct - n_det + 1:n_det, ] <- coef[["rate_det"]]
          }
          # The loadings come first among the coefficients.
          if (r > 0 & !is.null(coef[["rate_alpha"]])) {
            object[["priors"]][["a"]][["rate"]][1:n_alpha, ] <- coef[["rate_alpha"]]
          }
        }
      } else {
        v_i <- diag(coef[["v_i"]], tot_par)
        # Add priors for deterministic terms if they were specified
        if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
        }
      }
      object[["priors"]][["a"]][["v_inv"]] <- v_i
    }
    
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
  }
  
  ## Covar priors ----
  
  if (!structural & covar & k_endogen > 1) {
    
    n_covar <- k_endogen * (k_endogen - 1) / 2
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
      stop("Stochastic search variable selection is not implemented for error ",
           "correction sub-models. Consider using 'bvs' instead.")
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
  
  # Error term ----
  if (sv) {
    
    object <- .add_priors_sv_helper(object, sigma, k_endogen)
    
  } else {
    pos <- "u_sigma"
    if (error_prior == "wishart") {
      object[["priors"]][[pos]][["type"]] <- "wishart"
      help_df <- sigma[["df"]]
      object[["priors"]][[pos]][["df"]] <- NA_real_
      object[["priors"]][[pos]][["scale"]] = diag(sigma[["scale"]], k_endogen)
    }
    if (error_prior == "gamma") {
      object[["priors"]][[pos]][["type"]] <- "gamma"
      help_df <- sigma[["shape"]]
      object[["priors"]][[pos]][["shape"]] <- NA_real_
      object[["priors"]][[pos]][["rate"]] = matrix(sigma[["rate"]], k_endogen)
    }
    
    if (minnesota & !is.null(object[["data"]][["train"]][["x"]])) {
      # Store LS estimate of variance coviariance matrix for analytical solution
      object[["priors"]][[pos]][["u_sigma_inv"]] = minn[["sigma_inv"]]
    }
    
    if ("character" %in% class(help_df)) {
      if (grepl("k", help_df)) {
        # Transform character specification to expression and evaluate. The
        # number of endogenous variables is called k_endogen here, so 'k' is
        # supplied rather than looked up, where it does not exist.
        help_df <- eval(parse(text = help_df), envir = list(k = k_endogen))
      } else {
        stop("Use no other letter than 'k' in 'sigma$df' to indicate the number of endogenous variables.")
      }
    }
    
    if (help_df < 0) {
      stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
    }
    
    # Stored as given. BayesTS adds the rank to the posterior degrees of freedom
    # of the constant VECs with a Wishart prior itself, as the prior of the
    # loadings given the error covariance requires (Koop et al., 2010, eq. 8);
    # adding it here too counted it twice, and neither the gamma shape nor the
    # time varying models call for it. bvartools stores it as given as well.

    if (error_prior == "wishart") {
      object[["priors"]][[pos]][["df"]] <- help_df
    }
    if (error_prior == "gamma") {
      object[["priors"]][[pos]][["shape"]] <- matrix(help_df, k_endogen)
    }
  }
  
  return(object)
}