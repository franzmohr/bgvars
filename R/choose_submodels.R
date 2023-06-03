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
      pos <- which(names_object == names(criteria)[i])[which.min(criteria[[i]][, ic])]
      result[[i]] <- object[[pos]]
    }
    # Add rank selection here
    
  }
  names(result) <- names(criteria)
  
  class(result) <- c("bgvarest", "list")
  return(result)
}