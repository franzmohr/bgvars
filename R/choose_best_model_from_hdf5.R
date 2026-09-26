#' Choose Best Sub-Model
#'
#' Obtains model selection criteria for each model and gives
#' a table of the best performing model for each sub-model.
#'
#' @param folder path to a folder holding an export, that is, a \code{model.h5}
#' and a \code{submodels} directory beside it.
#' @param criterion selection criterion used to choose the best model using
#' \code{\link[bvartools]{selection_criteria}}.
#' @param plugin the point at which \code{"AIC"}, \code{"BIC"} and \code{"HQ"}
#' are evaluated, passed to \code{\link[bvartools]{selection_criteria}}.
#' \code{"mean"}, the default, uses the posterior mean and reads every model
#' file in full. \code{"best"} uses the draw with the highest log-likelihood,
#' which the light read already provides, so the point criteria then cost no
#' more than \code{"WAIC"} does. It is comparable only across models whose
#' chains are of the same length.
#'
#' @details
#' For the criteria that are calculated from the draws of the log-likelihood
#' alone (\code{"LL"}, \code{"WAIC"}, \code{"LOOIC"}) and for those calculated
#' from forecast errors, only the parts of a model file that enter the criterion
#' are read, that is, the draws of the log-likelihood, the forecast errors and
#' the handful of model specifications the criteria are calculated from. The
#' posterior draws of the coefficients make up most of an exported model and are
#' not touched, so the function does not need to hold a whole sub-model in
#' memory to rank it.
#'
#' \code{"AIC"}, \code{"BIC"} and \code{"HQ"} are evaluated at the posterior mean
#' of the parameters, see \code{\link[bvartools]{selection_criteria}}, which
#' needs the draws of the coefficients and the data. For them each model file is
#' read in full, one at a time, which takes longer. With
#' \code{plugin = "best"} they are evaluated at the draw with the highest
#' log-likelihood instead, which the light read already provides, and no model
#' file is read in full.
#'
#' A model file that holds neither draws of the log-likelihood nor forecast
#' errors -- a model an external sampler has not estimated, or failed on -- is
#' skipped with a warning that names it, and the positions of the remaining
#' models are those among all models of their sub-model. A sub-model none of
#' whose models has been estimated gets the position \code{NA}.
#'
#' @return A data frame containing the postion of the best performing model
#' for each sub-model.
#'
#' @export
choose_best_model_from_hdf5 <- function(folder, criterion, plugin = c("mean", "best")) {

  plugin <- match.arg(plugin)

  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " does not exist."))
  }

  filename_model <- file.path(folder, "model.h5")
  if (!file.exists(filename_model)) {
    stop(paste0("Folder ", folder, " does not contain a model.h5."))
  }

  model_file <- hdf5r::h5file(filename_model, mode = "r")
  on.exit(if (model_file$is_valid) model_file$close_all(), add = TRUE)

  # Sub-models ----
  #
  # Through the manifest rather than through the directory tree: which
  # sub-model a file belongs to is recorded there, and a folder that has been
  # moved or renamed since must not change the answer.
  manifest <- NULL
  if ("submodels" %in% names(model_file)) {
    manifest <- model_file[["submodels"]]$read()
  }
  model_file$close_all()

  if (is.null(manifest)) {
    stop("No submodels available.")
  }

  # Evaluated at the posterior mean, which the inputs read by
  # .read_selection_data_from_hdf5 cannot provide.
  # LML is not among them: it is read off the file, not evaluated at anything.
  #
  # With plugin = "best" the deviance is taken at the draw with the highest
  # log-likelihood, which the draws of the log-likelihood already hold, so even
  # AIC, BIC and HQ need nothing beyond what the light read provides and the
  # whole file is left alone.
  point_criterion <- criterion %in% c("AIC", "BIC", "HQ") && plugin == "mean"

  submodel <- unique(manifest[, "submodel"])

  result <- data.frame("submodel" = submodel, "position" = NA_integer_)
  for (i in seq_along(submodel)) {

    files <- manifest[manifest[, "submodel"] == submodel[i], "file"]

    # The criteria are obtained file by file, so that at no point more than one
    # model's draws are in memory. An element left NULL is a model without
    # results; choose_best_model() cannot choose it but keeps its position.
    sc <- vector(mode = "list", length = length(files))
    skipped <- character(0)
    for (j in seq_along(files)) {
      filename <- file.path(folder, files[j])
      object <- .read_selection_data_from_hdf5(filename)

      if (is.null(object[["posterior"]][["loglik"]]) &&
          is.null(object[["posterior"]][["forecast_errors"]])) {
        skipped <- c(skipped, files[j])
        next
      }

      if (point_criterion) {
        object <- bvartools::read_model_from_hdf5(filename = filename)
      }

      sc[[j]] <- .selection_criteria_from_hdf5(object, point_criterion, plugin)
    }

    if (length(skipped) > 0) {
      warning("Sub-model ", submodel[i], ": ", length(skipped), " of ", length(files),
              " model files hold neither draws of the log-likelihood nor forecast ",
              "errors and were skipped: ", paste0(skipped, collapse = ", "), ".",
              call. = FALSE)
    }

    if (length(skipped) == length(files)) {
      warning("Sub-model ", submodel[i], " has no estimated model, so its position is NA.",
              call. = FALSE)
      next
    }

    class(sc) <- c("selcritlist", "list")

    # Ties are broken by taking the first model.
    result[i, "position"] <- bvartools::choose_best_model(sc, criterion = criterion)[1]
  }

  return(result)
}

# Selection criteria of a model read from an export.
#
# A model read with .read_selection_data_from_hdf5 carries no coefficient
# draws, so bvartools cannot evaluate its log-likelihood at the posterior mean
# and warns that AIC, BIC and HQ are unavailable. Those criteria are not what
# such a model is read for -- a point criterion reads the whole file instead --
# so that warning, and only that one, is muffled. Once per file it would
# otherwise bury any warning that matters.
.selection_criteria_from_hdf5 <- function(object, point_criterion, plugin = "mean") {

  if (point_criterion) {
    return(bvartools::selection_criteria(object, plugin = plugin))
  }

  withCallingHandlers(
    bvartools::selection_criteria(object, plugin = plugin),
    warning = function(w) {
      if (grepl("could not be evaluated at the posterior mean", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
}

# Whether 'path' exists in an open HDF5 file.
#
# H5Lexists, which hdf5r's exists() calls, fails rather than answering FALSE
# when a group on the way to the last element is missing -- "posterior/loglik"
# of a model that has no posterior yet -- so the path is checked element by
# element.
.hdf5_path_exists <- function(h5, path) {

  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  for (i in seq_along(parts)) {
    if (!h5$exists(paste0(parts[seq_len(i)], collapse = "/"))) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Read the Inputs of a Selection Criterion from an HDF5 File
#'
#' Reads the parts of an exported model that
#' \code{\link[bvartools]{selection_criteria}} works on and returns them as a
#' model object of the class the file was written with.
#'
#' @param filename path to a single exported model.
#'
#' @details
#' The training data only enters the criteria through its dimensions, so its
#' shape is taken from the file and the data itself is not read. The same holds
#' for the specification of the model, of which only the entries a criterion is
#' calculated from are read.
#'
#' Handles are closed one by one instead of through \code{close_all}, which
#' walks every handle the file has ever opened and, over the hundreds of files
#' a global model consists of, costs more time than reading the draws.
#'
#' @return A list of the class recorded in the file.
#'
#' @noRd
.read_selection_data_from_hdf5 <- function(filename) {

  if (!file.exists(filename)) {
    stop(paste0("Model file ", filename, " does not exist."))
  }

  h5 <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5$is_valid) h5$close_all(), add = TRUE)
  handles <- list()

  # Specification ----
  #
  # Everything else the group records describes the model rather than entering
  # a criterion, and is available from the manifest for the models it is
  # wanted for.
  group <- h5[["model"]]
  handles <- c(handles, group)
  # 'algorithm' is here for the discounted models alone, whose criterion is the
  # exact log marginal likelihood their file carries rather than anything
  # estimated from a chain: selection_criteria() tells them apart by name, and
  # without it a discounted model would be scored as though it were a sampler.
  used_attrs <- c("k", "p", "m", "s", "n", "rank", "h", "structural", "rclass",
                  "algorithm")
  model <- list()
  for (attr_name in intersect(used_attrs, hdf5r::h5attr_names(group))) {
    model[[attr_name]] <- hdf5r::h5attr(group, attr_name)
  }

  object <- list("model" = model)

  # Training data ----
  #
  # Only the number of observations and the number of regressors are used.
  if (.hdf5_path_exists(h5, "data/train/y")) {
    dataset <- h5[["data/train/y"]]
    handles <- c(handles, dataset)
    object[["data"]][["train"]][["y"]] <- matrix(NA_real_, dataset$dims[1], 1L)
  }
  if (.hdf5_path_exists(h5, "data/train/x")) {
    dataset <- h5[["data/train/x"]]
    handles <- c(handles, dataset)
    object[["data"]][["train"]][["x"]] <- matrix(NA_real_, 1L, dataset$dims[2])
  }
  # The width of the error correction term, which is what the free parameters
  # of a cointegration matrix of a given rank are counted from.
  if (.hdf5_path_exists(h5, "data/train/w")) {
    dataset <- h5[["data/train/w"]]
    handles <- c(handles, dataset)
    object[["data"]][["train"]][["w"]] <- matrix(NA_real_, 1L, dataset$dims[2])
  }

  # Posterior ----
  if (.hdf5_path_exists(h5, "posterior/loglik")) {
    dataset <- h5[["posterior/loglik"]]
    handles <- c(handles, dataset)
    # Shaped from what the file says rather than from what hdf5r hands back: a
    # dataset with one row comes back as a plain vector, dimensions dropped, and
    # the criteria would then read one draw over the sample as one column of
    # many draws. A discounted model's log-likelihood is exactly one row.
    object[["posterior"]][["loglik"]] <-
      matrix(dataset$read(), nrow = dataset$dims[1], ncol = dataset$dims[2])
  }
  if (.hdf5_path_exists(h5, "posterior/forecast_errors")) {
    dataset <- h5[["posterior/forecast_errors"]]
    handles <- c(handles, dataset)
    object[["posterior"]][["forecast_errors"]] <- dataset$read()
    # Forecast errors are the only output that is labelled by variable.
    endogen <- hdf5r::h5attr(group, "endogen")
    object[["data"]][["original"]][["endogen"]] <-
      matrix(NA_real_, 1L, length(endogen), dimnames = list(NULL, endogen))
  }

  for (handle in handles) {
    handle$close()
  }
  h5$close()

  class(object) <- model[["rclass"]]

  return(object)
}
