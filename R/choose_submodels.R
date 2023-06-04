#' Susbmodel Selection
#' 
#' Selects the best country-specific submodel of a GVAR model.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{draw_posterior}}.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are "AIC", "BIC" (default) and "HQ".
#' @param select a character specifying how the best country model is selected.
#' If \code{"order"} (default), the country model with the overall minimum value of the
#' specified information criterion per country is selected. If \code{"rank"}, the
#' function selects the model, after which the selected information criterion increases
#' for the first time.
#' 
#' @export
choose_submodels <- function(object, ic = "BIC", select = "order") {
  
  obj_class <- class(object)
  
  if (!ic %in% c("AIC", "BIC", "HQ")) {
    stop("Invalid specification of argument 'ic'.")
  }
  
  # Obtain test statistics
  criteria <- submodel_test_satistics(object)
  names_object <- names(object)
  
  # Select final country models
  result <- NULL
  for (i in 1:length(criteria)) {
    # Order selection
    if (select == "order") {
      sub_pos <- which.min(criteria[[i]][, ic])
    }
    # Add rank selection here
    if (select == "rank") {
      crit <- criteria[[i]][, ic]
      sub_pos <- 1
      if (length(crit) > 1) {
        crit <- crit[-1] - crit[-length(crit)]
        crit <- which(crit < 0)
        if (length(crit) > 0) {
          sub_pos <- crit[1] + 1
        }
      }
      pos <- which(names_object == names(criteria)[i])[sub_pos]
      result[[i]] <- object[[pos]]
    }
  }
  names(result) <- names(criteria)
  
  class(result) <- obj_class
  return(result)
}