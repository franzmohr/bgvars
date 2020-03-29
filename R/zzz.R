# Load needed functions from different packages for fast access
#' @importFrom bvartools bvs
# @importFrom bvartools kalman_dk
#' @importFrom bvartools loglik_normal
#' @importFrom bvartools post_coint_kls
#' @importFrom bvartools post_coint_kls_sur
#' @importFrom bvartools post_normal
#' @importFrom bvartools post_normal_sur
#' @importFrom bvartools ssvs
# @importFrom stats rgamma
#' @importFrom stats rWishart

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("bgvars", libpath)
}