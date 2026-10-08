# Shared implementation behind write_to_hdf5.gvarmodel and write_to_hdf5.gvecmodel.
#
# Layout produced:
#
#   <folder>/model.h5              global data, weight matrices and the manifest
#   <folder>/submodels/<name>/001.h5, 002.h5, ...
#
# Everything written once by one process goes into model.h5. Each sub-model
# keeps its own file, because bayests::write_to_hdf5 writes a model at the
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

# Attaches 'value' to 'object' as the attribute 'name', and closes it again.
#
# hdf5r's h5attr<- creates the attribute and leaves its handle open for the
# garbage collector. HDF5 keeps a file open for as long as anything in it is,
# and on Windows the file is then locked against every other process, BayesTS
# among them. The export ends in close_all(), which releases such handles, so
# this is not what keeps a file open today; closing each attribute as it is
# written means that stays true if the file comes to be closed another way.
# bvartools writes its attributes the same way. An attribute that is already
# there is replaced.
.write_attr <- function(object, name, value) {
  if (object$attr_exists(name)) {
    object$attr_delete(name)
  }

  attribute <- object$create_attr(name, robj = value)
  attribute$close()

  invisible(NULL)
}

# Writes a time-series object together with the two attributes that make it one
# again on the way back: its variable names and its tsp.
.write_series <- function(group, name, x) {
  group[[name]] <- x
  dataset <- group[[name]]
  on.exit(dataset$close(), add = TRUE)
  .write_attr(dataset, "variables", dimnames(x)[[2]])
  .write_attr(dataset, "tsp", stats::tsp(x))

  invisible(NULL)
}

# One row per sub-model file, describing the model in it.
#
# This is what the file names used to carry. Everything here is read back off
# the model's own specification, so the manifest cannot disagree with the file
# it points at.
# One row of the manifest: what a model's specification says about it, beside
# the sub-model it belongs to and the file it is written to. The same row is
# built when a model is written and when map_submodels() has changed one on
# disk, so the two cannot describe a model differently.
.manifest_row <- function(specs, submodel, file) {

  field <- function(specs, name, mode) {
    value <- specs[[name]]
    if (is.null(value) || length(value) != 1) {
      return(as.vector(NA, mode = mode))
    }
    as.vector(value, mode = mode)
  }

  data.frame(
    submodel = submodel,
    file = file,
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
    # A model that keeps every draw carries no 'thin', so that is one rather
    # than missing: NA would read like an export from before the manifest
    # recorded thinning.
    thin = if (is.null(specs[["thin"]])) 1L else field(specs, "thin", "integer"),
    # The discounted models only. One is the value at which the quantity a
    # discount governs does not move, so a sub-model that carries neither is
    # recorded as one rather than as missing -- which is what every other
    # algorithm is: coefficients and an error covariance that are what they are.
    delta_beta = if (is.null(specs[["delta_beta"]])) 1 else field(specs, "delta_beta", "numeric"),
    delta_sigma = if (is.null(specs[["delta_sigma"]])) 1 else field(specs, "delta_sigma", "numeric"),
    stringsAsFactors = FALSE
  )
}

.submodel_manifest <- function(object) {

  rows <- list()

  for (submodel in names(object[["submodels"]])) {
    models <- object[["submodels"]][[submodel]]

    for (i in seq_along(models)) {
      rows[[length(rows) + 1]] <- .manifest_row(
        models[[i]][["model"]], submodel,
        .submodel_file(i = i, submodel = submodel, n = length(models)))
    }
  }

  if (length(rows) == 0) {
    return(NULL)
  }

  result <- do.call("rbind", rows)
  rownames(result) <- NULL

  return(result)
}

# The manifest of an append: the rows on disk, without the sub-models the batch
# brings, followed by the rows of the batch.
#
# A manifest read back from HDF5 comes as a data frame of the same columns, but
# hdf5r does not promise the storage types it read them into, so the kept rows
# are coerced to the types of the batch before the two are stacked.
.merged_manifest <- function(model_file, manifest) {

  if (!model_file$exists("submodels")) {
    return(manifest)
  }

  existing <- as.data.frame(model_file[["submodels"]]$read(), stringsAsFactors = FALSE)
  if (is.null(manifest)) {
    return(existing)
  }

  missing <- setdiff(names(manifest), names(existing))
  if (length(missing) > 0) {
    stop("The export in this folder has no column ", paste(missing, collapse = ", "),
         " in its manifest, so it was written by another version of the package ",
         "and cannot be appended to.")
  }

  kept <- existing[!existing[["submodel"]] %in% manifest[["submodel"]],
                   names(manifest), drop = FALSE]
  for (column in names(manifest)) {
    mode(kept[[column]]) <- mode(manifest[[column]])
  }

  result <- rbind(kept, manifest)
  rownames(result) <- NULL

  result
}

# Path of one sub-model's file, relative to the export folder.
#
# 'n' is the number of models of the sub-model. The numbers are padded to the
# width of the largest, and to three digits at least, so that the files sort as
# text in the order they were written: a sub-model of a thousand models would
# otherwise put 1000.h5 before 101.h5 in any directory listing, and with it
# bayests::read_models_from_folder.
.submodel_file <- function(i, submodel, n = i) {
  width <- max(3L, nchar(as.character(n)))
  paste0("submodels/", submodel, "/", formatC(i, width = width, flag = "0"), ".h5")
}

# How many cores an export uses when the caller does not say.
#
# Capped rather than taken whole: an export is not the only thing the machine
# is doing, and the gain flattens out well before the last core. Eight was
# where it flattened on a sixteen-core machine.
.default_cores <- function(limit = 8L) {

  # R CMD check sets this, and CRAN allows a check no more than two cores.
  # Honoured here
  # rather than in every example and test, so that no call site has to know.
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    limit <- 2L
  }

  available <- parallel::detectCores()

  # detectCores() is documented to return NA where it cannot tell.
  if (is.na(available)) {
    return(1L)
  }

  as.integer(max(1L, min(limit, available)))
}

# Writes every model of one sub-model.
#
# The unit a worker is given when an export runs on more than one core, and the
# body of the loop when it does not. It takes the file names rather than
# working them out, so that a worker is sent nothing but the models it writes
# and the paths to write them to.
.write_submodel_to_hdf5 <- function(task) {

  # Through bvartools, which re-exports the generic of bayests: a worker loads
  # the namespace a call names, and the methods for the sub-models' VAR and
  # VEC classes are registered by bvartools. bayests:: would load the generic
  # without them.
  for (i in seq_along(task[["models"]])) {
    bvartools::write_to_hdf5(task[["models"]][[i]],
                             filename = task[["filenames"]][i])
  }

  invisible(NULL)
}

.write_gvar_to_hdf5 <- function(object, folder, overwrite = FALSE,
                                mc.cores = .default_cores(), append = FALSE) {

  if (!dir.exists(folder)) {
    stop(paste0("Folder ", folder, " does not exist."))
  }

  filename_model <- file.path(folder, "model.h5")
  path_submodels <- file.path(folder, "submodels")

  # Checked before anything is written, so that the common mistake -- exporting
  # twice into the same folder -- costs nothing and leaves the first export
  # intact. bayests::write_to_hdf5 refuses an existing file as well, and
  # finding that out on sub-model 57 of 132 is what the check up here avoids.
  if (append) {
    if (!file.exists(filename_model)) {
      stop(paste0(filename_model, " does not exist, so there is no export to append to."))
    }
  } else {
    for (path in c(filename_model, path_submodels)) {
      if (file.exists(path)) {
        if (!overwrite) {
          stop(paste0(path, " already exists. Use 'overwrite = TRUE' to replace it."))
        }
        unlink(path, recursive = TRUE)
      }
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

  manifest <- .submodel_manifest(object)

  if (append) {

    # The global data, the index and the weight matrices describe the global
    # model rather than the batch, and every batch of one model agrees on them,
    # so an append leaves them as the first write left them and touches the
    # manifest alone.
    manifest <- .merged_manifest(model_file, manifest)
    if (!is.null(manifest)) {
      if (model_file$exists("submodels")) {
        model_file$link_delete("submodels")
      }
      model_file[["submodels"]] <- manifest
    }

  } else {

    .write_attr(model_file, "rclass", class(object))

    group_global <- model_file$create_group("global")
    for (i in c("endogen", "exogen", "deterministic")) {
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

    if (!is.null(manifest)) {
      model_file[["submodels"]] <- manifest
    }
  }

  model_file$close_all()

  # Sub-models ----
  #
  # One file per model, named by position rather than by specification: the
  # manifest above says which model each file holds.
  if (!is.null(manifest)) {
    dir.create(path_submodels, showWarnings = FALSE)
    for (submodel in names(object[["submodels"]])) {
      # A sub-model the export already holds is replaced as a whole: its files
      # are numbered by position, and a batch with fewer models than the one on
      # disk would otherwise leave the surplus behind for the manifest to
      # point at.
      path <- file.path(path_submodels, submodel)
      if (dir.exists(path)) {
        unlink(path, recursive = TRUE)
      }
      dir.create(path)
    }

    # One task per sub-model, carrying only the models it writes and where they
    # go, so that a worker is never sent the whole object.
    tasks <- lapply(names(object[["submodels"]]), function(submodel) {
      models <- object[["submodels"]][[submodel]]
      list("models" = models,
           "filenames" = file.path(folder,
                                   vapply(seq_along(models), .submodel_file,
                                          character(1), submodel = submodel,
                                          n = length(models))))
    })

    # No more workers than there is work for them, and never a cluster to write
    # a single sub-model: starting one costs more than it can save there.
    if (is.null(mc.cores) || is.na(mc.cores)) {
      mc.cores <- 1L
    }
    mc.cores <- min(as.integer(mc.cores), length(tasks))

    if (mc.cores < 2L) {
      lapply(tasks, .write_submodel_to_hdf5)
    } else {
      # A socket cluster rather than parallel::mclapply, which is what the rest
      # of the package uses: mclapply forks, has nothing to fall back on under
      # Windows and runs serially there, and an export is long enough that
      # quietly doing nothing is the wrong trade. Every sub-model has its own
      # files and its directory already exists, so the workers share nothing.
      cluster <- parallel::makePSOCKcluster(mc.cores)
      on.exit(parallel::stopCluster(cluster), add = TRUE)

      # Balanced, because sub-models differ in how many models they hold and in
      # how large those are.
      parallel::parLapplyLB(cluster, tasks, .write_submodel_to_hdf5)
    }
  }

  invisible(folder)
}
