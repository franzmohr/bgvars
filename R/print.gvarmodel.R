#' Printing Model Information
#'
#' print methods for objects of class 'gvarmodel' and 'gvecmodel'.
#'
#' @param x an object of class 'gvarmodel' or 'gvecmodel'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Printed are the variables of every sub-model, with a mark for each variable
#' that the sub-model contains, the global variables, whether weight matrices
#' have been added and, once sub-models have been added, a line per sub-model:
#' how many candidate models it holds, whether they have been estimated, and
#' the specification of the model when there is only one.
#'
#' @return The object that was printed, invisibly.
#'
#' @export
print.gvarmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  .print_global_model(x, "Global Vector Autoregressive Model", digits = digits, ...)
}

#' @rdname print.gvarmodel
#' @export
print.gvecmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  .print_global_model(x, "Global Vector Error Correction Model", digits = digits, ...)
}

# Shared by both print methods. The method for 'gvarmodel' stopped with a note
# to update it as soon as a model held sub-models, which is when a model is
# printed most, and a 'gvecmodel' had no method and was printed as a raw list.
.print_global_model <- function(x, title, digits, ...) {

  # Used variables
  avail_vars <- x[["global"]][["index"]][, c("submodel", "variable")]

  submodels <- unique(avail_vars[, "submodel"])
  vars <- unique(avail_vars[, "variable"])

  availability_matrix <- matrix("", nrow = length(submodels), ncol = length(vars))
  dimnames(availability_matrix) <- list(submodels, vars)

  for (i in submodels) {
    vars_i <- avail_vars[avail_vars[, "submodel"] == i, "variable"]
    availability_matrix[i, vars_i] <- "x"
  }

  res <- data.frame(submodel = submodels, as.data.frame(availability_matrix, row.names = FALSE))

  cat("\n", title, "\n", sep = "")

  cat("\nVariables in sub-models:\n\n")
  print(res, digits = digits, row.names = FALSE, ...)

  if (!is.null(x[["global"]][["exogen"]])) {
    res <- dimnames(x[["global"]][["exogen"]])[[2]]
    cat("\nGlobal variables:\n\n")
    cat("\t", paste0(res, collapse = "\n\t"), "\n", sep = "")
  }

  cat("\nWeight matrices: ",
      if (is.null(x[["weights"]])) "not added" else "added", "\n", sep = "")

  if (length(x[["submodels"]]) > 0) {

    rows <- lapply(names(x[["submodels"]]), function(i) {
      models <- x[["submodels"]][[i]]
      estimated <- vapply(models, function(m) !is.null(m[["posterior"]]), logical(1))
      row <- data.frame(submodel = i,
                        models = length(models),
                        estimated = if (all(estimated)) "yes" else if (any(estimated)) "partly" else "no",
                        stringsAsFactors = FALSE)
      if (length(models) == 1) {
        specs <- models[[1]][["model"]]
        row[["type"]] <- specs[["type"]]
        row[["p_endogen"]] <- specs[["p_endogen"]]
        row[["p_exogen"]] <- specs[["p_exogen"]]
        row[["rank"]] <- if (is.null(specs[["rank"]])) NA_integer_ else specs[["rank"]]
      } else {
        row[["type"]] <- NA_character_
        row[["p_endogen"]] <- NA_integer_
        row[["p_exogen"]] <- NA_integer_
        row[["rank"]] <- NA_integer_
      }
      row
    })
    res <- do.call(rbind, rows)
    # Level models have no cointegration rank.
    if (all(is.na(res[["rank"]]))) {
      res[["rank"]] <- NULL
    }

    cat("\nSub-models:\n\n")
    print(res, digits = digits, row.names = FALSE, ...)
  }

  invisible(x)
}
