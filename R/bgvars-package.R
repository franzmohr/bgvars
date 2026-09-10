#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib bgvars, .registration = TRUE
#' @importFrom bvartools add_initial_values
#' @importFrom bvartools add_posterior_coefficients
#' @importFrom bvartools add_posterior_forecasts
#' @importFrom bvartools add_posterior_loglik
#' @importFrom bvartools add_priors
#' @importFrom bvartools align_model_obs
#' @importFrom bvartools choose_best_model
#' @importFrom bvartools get_model_specifications
#' @importFrom bvartools inclusion_prior
#' @importFrom bvartools irf
#' @importFrom bvartools minnesota_prior
#' @importFrom bvartools rescale_error_correction
#' @importFrom bvartools scale_error_correction
#' @importFrom bvartools selection_criteria
#' @importFrom bvartools ssvs_prior
#' @importFrom bvartools use_expanding_window
#' @importFrom bvartools vec_to_var
#' @importFrom bvartools write_to_hdf5
#' @importFrom coda thin
#' @importFrom Rcpp sourceCpp
#' @importFrom stats predict
#' @importFrom stats window
#' @exportPattern "^[[:alpha:]]+"
## usethis namespace: end
NULL
