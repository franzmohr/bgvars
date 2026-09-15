#' Thinning Posterior Draws
#'
#' Forwards model input to the same function for individual models.
#'
#' @param x an object of class 'gvarmodel'.
#' @param thin an integer specifying the thinning interval between successive
#' values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Every model of every sub-model is thinned by the method of
#' \pkg{bvartools}, which counts the labels of the draws from the ones they
#' already carry. A sub-model whose sampler kept one draw in \code{k} is
#' therefore labelled as keeping one in \code{k * thin} afterwards. Every
#' sub-model has to be estimated.
#'
#' @return An object of class 'gvarmodel'.
#'
#' @export
#' @method thin gvarmodel
thin.gvarmodel <- function(x, thin = 10, ...) {
  x[["submodels"]] <- .thin_submodels(x[["submodels"]], thin = thin, ...)
  return(x)
}


#' Thinning Posterior Draws
#'
#' Forwards model input to the same function for individual models.
#'
#' @param x an object of class 'gvecmodel'.
#' @param thin an integer specifying the thinning interval between successive
#' values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Every model of every sub-model is thinned by the method of
#' \pkg{bvartools}, which counts the labels of the draws from the ones they
#' already carry. A sub-model whose sampler kept one draw in \code{k} is
#' therefore labelled as keeping one in \code{k * thin} afterwards. Every
#' sub-model has to be estimated.
#'
#' @return An object of class 'gvecmodel'.
#'
#' @export
#' @method thin gvecmodel
thin.gvecmodel <- function(x, thin = 10, ...) {
  x[["submodels"]] <- .thin_submodels(x[["submodels"]], thin = thin, ...)
  return(x)
}


# Thins every model of every sub-model, keeping the class of the list each
# sub-model's models are held in.
.thin_submodels <- function(submodels, thin, ...) {
  for (i in names(submodels)) {
    for (j in seq_along(submodels[[i]])) {
      if (is.null(submodels[[i]][[j]][["posterior"]])) {
        stop("Sub-model ", i, " has not been estimated, so it has no draws to thin.")
      }
      submodels[[i]][[j]] <- coda::thin(submodels[[i]][[j]], thin = thin, ...)
    }
  }
  submodels
}
