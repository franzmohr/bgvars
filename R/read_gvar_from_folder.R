#' Import a Global Model from a Folder
#'
#' Reads back a global model that was exported with \code{\link[bvartools]{write_to_hdf5}},
#' including any posterior draws an external sampler has written into the
#' sub-model files since.
#'
#' @param folder path to a folder holding an export, that is, a \code{model.h5}
#' and a \code{submodels} directory beside it.
#'
#' @details
#'
#' The sub-models are found through the \code{submodels} manifest in
#' \code{model.h5}, which records the file each of them was written to. The
#' directory tree is not consulted: a name a sub-model is known by --- which
#' country or region it is --- is information the writer recorded, and reading
#' it back off a path would make it depend on how the folder was moved or
#' renamed since.
#'
#' The class of the result is likewise the one the export recorded, so a
#' 'gvecmodel' comes back as one.
#'
#' @return An object of class 'gvarmodel' or 'gvecmodel'.
#'
#' @examples
#'
#' # Load data
#' data("gvar2023")
#' submodel_data <- gvar2023[["submodel_data"]]
#' global_data <- gvar2023[["global_data"]]
#' 
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#'
#' # Set up a model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#'                            
#' # Generate and add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 2014:2016)
#'                               
#' # Create sub-models
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp", "r"),
#'                         p_endogen = 1,
#'                         exogen = c("y", "Dp", "r"),
#'                         p_exogen = 1,
#'                         global = "poil",
#'                         s = 1,
#'                         error = "wishart",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 0.0001))
#'                      
#' # Add initial values
#' object <- add_initial_values(object)
#'
#' folder <- file.path(tempdir(), "bgvars-example-read-gvar")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE)
#'
#' # Read it back
#' object <- read_gvar_from_folder(folder)
#'
#' @export
read_gvar_from_folder <- function(folder, submodels = NULL) {

  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " does not exist."))
  }

  filename_model <- file.path(folder, "model.h5")
  if (!file.exists(filename_model)) {
    stop(paste0("Folder ", folder, " does not contain a model.h5."))
  }

  model_file <- hdf5r::h5file(filename_model, mode = "r")
  on.exit(if (model_file$is_valid) model_file$close_all(), add = TRUE)

  result <- list("global" = list(),
                 "weights" = NULL,
                 "submodels" = NULL)

  # Global data ----
  if (!"global" %in% names(model_file)) {
    stop(paste0(filename_model, " does not contain global data."))
  }

  group_global <- model_file[["global"]]
  for (i in c("endogen", "exogen", "deterministic")) {
    if (i %in% names(group_global)) {
      result[["global"]][[i]] <- .read_series(group_global[[i]])
    }
  }
  result[["global"]][["index"]] <- group_global[["index"]]$read()

  # Weights ----
  #
  # In the order of the sub-models in the index, which is the order
  # add_weight_matrices() puts them in. HDF5 hands its group members back
  # alphabetically, so taking them as they come would quietly reorder the list
  # against the model it belongs to.
  if ("weights" %in% names(model_file)) {
    group_weights <- model_file[["weights"]]

    stored <- names(group_weights)
    ordered <- unique(result[["global"]][["index"]][, "submodel"])
    ordered <- c(ordered[ordered %in% stored], setdiff(stored, ordered))

    result[["weights"]] <- lapply(ordered, function(i) group_weights[[i]]$read())
    names(result[["weights"]]) <- ordered
  }

  # Sub-models ----
  #
  # Through the manifest rather than through the directory tree: which
  # sub-model a file belongs to is recorded there, and a folder that has been
  # moved or renamed since must not change the answer.
  manifest <- NULL
  if ("submodels" %in% names(model_file)) {
    manifest <- model_file[["submodels"]]$read()
  }

  result_class <- hdf5r::h5attr(model_file, "rclass")
  model_file$close_all()

  if (!is.null(manifest) && nrow(manifest) > 0) {
    
    if (!is.null(submodels)) {
      manifest <- .select_manifest_rows(manifest, submodels)
    }
    
    result[["submodels"]] <- .read_submodels(folder, manifest)
  }

  class(result) <- result_class

  return(result)
}

# Restores a time-series object from a dataset written by .write_series.
.read_series <- function(dataset) {

  result <- stats::ts(as.matrix(hdf5r::readDataSet(dataset)),
                      class = c("mts", "ts", "matrix"))
  dimnames(result) <- list(NULL, hdf5r::h5attr(dataset, "variables"))
  stats::tsp(result) <- hdf5r::h5attr(dataset, "tsp")

  return(result)
}

# One 'modellist' per sub-model, in the order the manifest lists them.
.read_submodels <- function(folder, manifest) {

  # A group column is not written today -- every sub-model file holds one model
  # at its root -- but reading it when it is there keeps this working if that
  # changes.
  groups <- if ("group" %in% names(manifest)) {
    manifest[, "group"]
  } else {
    rep("", nrow(manifest))
  }
  groups[is.na(groups)] <- ""

  result <- list()

  for (submodel in unique(manifest[, "submodel"])) {

    rows <- which(manifest[, "submodel"] == submodel)

    models <- list()
    for (i in seq_along(rows)) {
      filename <- file.path(folder, manifest[rows[i], "file"])

      if (!file.exists(filename)) {
        stop(paste0("The manifest of this export names ", manifest[rows[i], "file"],
                    ", which is not in ", folder, "."))
      }

      # The group is only passed when there is one, so that an export in the
      # layout written here -- one model at the root of its file -- reads back
      # with a bvartools that predates the group argument.
      models[[i]] <- if (groups[rows[i]] == "") {
        bvartools::read_model_from_hdf5(filename = filename)
      } else {
        bvartools::read_model_from_hdf5(filename = filename,
                                        group = groups[rows[i]])
      }
    }

    class(models) <- c("modellist", "list")
    result[[submodel]] <- models
  }

  return(result)
}

# The rows of the manifest that a selection of one model per sub-model picks
# out.
#
# The 'position' of a model is its position among the models of its sub-model,
# which is the order the manifest lists them in, the order they were written
# in, and hence the order bvartools::read_models_from_folder reads them back
# in. Selections made from either route therefore mean the same thing.
.select_manifest_rows <- function(manifest, submodels) {

  if (!is.data.frame(submodels)) {
    submodels <- as.data.frame(submodels)
  }

  missing_cols <- setdiff(c("submodel", "position"), names(submodels))
  if (length(missing_cols) > 0) {
    stop("Argument 'submodels' must have the column(s) ",
         paste0(missing_cols, collapse = ", "), ".")
  }

  unknown <- setdiff(submodels[, "submodel"], manifest[, "submodel"])
  if (length(unknown) > 0) {
    stop("This export does not contain the sub-model(s) ",
         paste0(unknown, collapse = ", "), ".")
  }

  rows <- integer(0)
  for (i in seq_len(nrow(submodels))) {

    submodel <- submodels[i, "submodel"]
    position <- submodels[i, "position"]

    available <- which(manifest[, "submodel"] == submodel)

    if (is.na(position) || position < 1 || position > length(available)) {
      stop("Sub-model ", submodel, " was exported with ", length(available),
           " model(s), so position ", position, " does not exist.")
    }

    rows <- c(rows, available[position])
  }

  return(manifest[rows, , drop = FALSE])
}
