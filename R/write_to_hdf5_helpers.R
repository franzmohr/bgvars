# Shared implementation behind write_to_hdf5.gvarmodel and write_to_hdf5.gvecmodel.
#
# Layout produced:
#
#   <folder>/model.h5              global data, weight matrices and the manifest
#   <folder>/submodels/<name>/001.h5, 002.h5, ...
#
# Everything written once by one process goes into model.h5. Each sub-model
# keeps its own file, because bvartools::write_to_hdf5 writes a model at the
# root of a file it creates and takes no group, and because one file per model
# is what lets an external driver run several of them at once -- HDF5 has no
# support for concurrent writers to one file.
#
# The file names carry no model specification. The manifest in model.h5 does,
# so a caller can find the model it wants without opening and parsing a hundred
# file names, and two models whose specifications happen to agree cannot
# collide.

# Rows of a chunk that keep one chunk under `max_bytes` uncompressed.
#
# Chunk shape decides how well the weight matrices compress: they are mostly
# zeros in long runs along a row, so a chunk spanning the full width squeezes
# far better than the small square chunks hdf5r picks on its own. Whole rows,
# then, and as many of them as fit in a chunk worth holding in memory.
.hdf5_chunk_dims <- function(x, max_bytes = 8e6) {
  n_rows <- nrow(x)
  n_cols <- ncol(x)

  if (n_rows == 0 || n_cols == 0) {
    return(NULL)
  }

  rows_per_chunk <- max(1, floor(max_bytes / (8 * n_cols)))

  c(min(n_rows, rows_per_chunk), n_cols)
}

# Writes a numeric matrix as a compressed dataset.
.write_compressed <- function(group, name, x, gzip_level = 9) {
  chunk_dims <- .hdf5_chunk_dims(x)

  if (is.null(chunk_dims)) {
    group[[name]] <- x
    return(invisible(NULL))
  }

  group$create_dataset(name, x, chunk_dims = chunk_dims, gzip_level = gzip_level)

  invisible(NULL)
}

# Writes a time-series object together with the two attributes that make it one
# again on the way back: its variable names and its tsp.
.write_series <- function(group, name, x) {
  group[[name]] <- x
  hdf5r::h5attr(group[[name]], "variables") <- dimnames(x)[[2]]
  hdf5r::h5attr(group[[name]], "tsp") <- stats::tsp(x)

  invisible(NULL)
}

# One row per sub-model file, describing the model in it.
#
# This is what the file names used to carry. Everything here is read back off
# the model's own specification, so the manifest cannot disagree with the file
# it points at.
.submodel_manifest <- function(object) {
  field <- function(specs, name, mode) {
    value <- specs[[name]]
    if (is.null(value) || length(value) != 1) {
      return(as.vector(NA, mode = mode))
    }
    as.vector(value, mode = mode)
  }

  rows <- list()

  for (submodel in names(object[["submodels"]])) {
    models <- object[["submodels"]][[submodel]]

    for (i in seq_along(models)) {
      specs <- models[[i]][["model"]]

      rows[[length(rows) + 1]] <- data.frame(
        submodel = submodel,
        file = .submodel_file(submodel, i),
        algorithm = field(specs, "algorithm", "character"),
        type = field(specs, "type", "character"),
        k = field(specs, "k", "integer"),
        p_endogen = field(specs, "p_endogen", "integer"),
        k_exogen = field(specs, "k_exogen", "integer"),
        p_exogen = field(specs, "p_exogen", "integer"),
        m_global = field(specs, "m_global", "integer"),
        s_global = field(specs, "s_global", "integer"),
        n = field(specs, "n", "integer"),
        # VEC models only. NA for a VARX, which has no cointegration rank.
        rank = field(specs, "rank", "integer"),
        varsel = field(specs, "varsel", "character"),
        structural = field(specs, "structural", "logical"),
        tvp = field(specs, "tvp", "logical"),
        error = field(specs, "error", "character"),
        iterations = field(specs, "iterations", "integer"),
        burnin = field(specs, "burnin", "integer"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(NULL)
  }

  result <- do.call("rbind", rows)
  rownames(result) <- NULL

  return(result)
}

# Path of one sub-model's file, relative to the export folder.
.submodel_file <- function(submodel, i) {
  paste0("submodels/", submodel, "/", sprintf("%03d", i), ".h5")
}

.write_gvar_to_hdf5 <- function(object, folder, overwrite = FALSE) {

  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " does not exist."))
  }

  filename_model <- file.path(folder, "model.h5")
  path_submodels <- file.path(folder, "submodels")

  # Checked before anything is written, so that the common mistake -- exporting
  # twice into the same folder -- costs nothing and leaves the first export
  # intact. bvartools::write_to_hdf5 refuses an existing file as well, and
  # finding that out on sub-model 57 of 132 is what the check up here avoids.
  for (path in c(filename_model, path_submodels)) {
    if (file.exists(path)) {
      if (!overwrite) {
        stop(paste0(path, " already exists. Use 'overwrite = TRUE' to replace it."))
      }
      unlink(path, recursive = TRUE)
    }
  }

  # Global data and weights ----
  #
  # No try() around any of this. A failure here used to be swallowed and the
  # function returned as though it had worked, leaving a folder that looked
  # complete and was not.
  model_file <- hdf5r::h5file(filename_model, mode = "a")
  # Guarded, because the handle is closed explicitly below once the file is
  # complete and closing it twice is an error in hdf5r. This is here for the
  # paths that do not reach that point.
  on.exit(if (model_file$is_valid) model_file$close_all(), add = TRUE)

  hdf5r::h5attr(model_file, "rclass") <- class(object)

  group_global <- model_file$create_group("global")
  for (i in c("endogen", "exogen")) {
    if (!is.null(object[["global"]][[i]])) {
      .write_series(group_global, i, object[["global"]][[i]])
    }
  }
  group_global[["index"]] <- object[["global"]][["index"]]

  if (!is.null(object[["weights"]])) {
    group_weights <- model_file$create_group("weights")
    for (i in names(object[["weights"]])) {
      .write_compressed(group_weights, i, object[["weights"]][[i]])
    }
  }

  manifest <- .submodel_manifest(object)
  if (!is.null(manifest)) {
    model_file[["submodels"]] <- manifest
  }

  model_file$close_all()

  # Sub-models ----
  #
  # One file per model, named by position rather than by specification: the
  # manifest above says which model each file holds.
  if (!is.null(manifest)) {
    dir.create(path_submodels)

    for (submodel in names(object[["submodels"]])) {
      dir.create(file.path(path_submodels, submodel))

      models <- object[["submodels"]][[submodel]]
      for (i in seq_along(models)) {
        bvartools::write_to_hdf5(models[[i]],
                                 filename = file.path(folder, .submodel_file(submodel, i)))
      }
    }
  }

  invisible(folder)
}
