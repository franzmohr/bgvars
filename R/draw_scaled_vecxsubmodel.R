# Posterior simulation of one VECX sub-model on a scaled error correction term.
#
# The series in the error correction term of a GVEC sub-model are on the scale
# of the data, and those scales differ a great deal: a log level and an
# inflation rate of the same unit sit in the same term, and so does a linear
# trend running into the hundreds. The prior on the cointegration space is
# isotropic -- add_priors() builds it as diag(p_tau_i, k_ect) -- so it treats
# every one of those series as if they were comparably scaled, which they are
# not.
#
# scale_error_correction() divides each stochastic series by the standard
# deviation of its differences, and the trend by its own standard deviation,
# which puts them on a comparable footing for the sampler.
# rescale_error_correction() puts the draws back on the scale of the data
# afterwards. Both come from bvartools, through the 'bvecmodel' class that a
# 'vecxsubmodel' inherits.
#
# The transformation is exact: with D the diagonal matrix of scaling factors,
# the error correction term alpha beta' w is unchanged by writing it as
# alpha (D beta)' (D^-1 w), so what is estimated is the same model. Only beta
# is affected, alpha is not.
.draw_scaled_vecxsubmodel <- function(object, ...) {

  object <- scale_error_correction(object)

  factors <- attr(object[["data"]][["train"]][["w"]], "scale")

  # The starting value of beta was computed from the data on their own scale,
  # by add_initial_values(), and has to follow the same transformation as the
  # series it multiplies. It is put back afterwards: those estimates are what a
  # replication compares with published cointegrating vectors, so the object
  # must not be left holding a transformed copy of them.
  # The starting values are stored stacked, as a single column of k_ect * r
  # elements, so they are reshaped before the factors are applied and stacked
  # again afterwards. diag() is given an explicit size, because for a single
  # factor it would otherwise read its argument as the size of an identity
  # matrix rather than as the diagonal to build.
  initial_beta <- object[["initial"]][["beta"]]
  if (!is.null(initial_beta)) {
    scaling <- diag(factors, nrow = length(factors))
    object[["initial"]][["beta"]] <-
      matrix(scaling %*% matrix(initial_beta, length(factors)), ncol = 1)
  }

  object <- add_posterior_coefficients(object, ...)

  object <- rescale_error_correction(object)
  object[["initial"]][["beta"]] <- initial_beta

  return(object)
}
