# Estimating a global model through its export, one model at a time.
#
# What the methods for a 'gvarexport' have to deliver is the posterior a global
# model in memory delivers: the same priors, the same starting values, the same
# seeds and so the same draws, with nothing but one model per worker in memory
# while they run.

export_model <- function(type = "gvar") {
  object <- if (type == "gvar") {
    add_submodels(gvar_object(), endogen = c("y", "Dp"), p_endogen = 1,
                  exogen = c("y", "Dp"), p_exogen = 1, deterministic = "const",
                  iterations = 20, burnin = 5)
  } else {
    add_submodels(gvec_object(), endogen = c("y", "Dp"), p_endogen = 1,
                  exogen = c("y", "Dp"), p_exogen = 2, const = "unrestricted",
                  trend = "restricted", r = 0:1, iterations = 20, burnin = 5)
  }
  object
}

export_priors <- function(object, type) {
  if (type == "gvar") {
    add_priors(object, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))
  } else {
    add_priors(object, coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = 3, scale = 1e-4))
  }
}

written_export <- function(object, name) {
  folder <- file.path(tempdir(), paste0("bgvars-export-", name))
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
  read_export(folder)
}

draws_of <- function(x, submodel, position) {
  x[["submodels"]][[submodel]][[position]][["posterior"]][["a"]][["coeffs"]]
}

test_that("an export is estimated to the posterior of the model in memory", {
  for (type in c("gvar", "gvec")) {

    object <- export_model(type)
    export <- written_export(object, type)

    in_memory <- export_priors(object, type)
    in_memory <- add_seed(add_initial_values(in_memory), 1234)
    in_memory <- add_posterior_loglik(add_posterior_coefficients(in_memory))

    on_disk <- export_priors(export, type)
    on_disk <- add_initial_values(on_disk)
    on_disk <- add_seed(on_disk, 1234)
    on_disk <- add_posterior_coefficients(on_disk)
    on_disk <- add_posterior_loglik(on_disk)
    expect_s3_class(on_disk, "gvarexport")

    back <- read_gvar_from_folder(on_disk[["folder"]])
    expect_identical(names(back[["submodels"]]), names(in_memory[["submodels"]]))
    for (submodel in names(in_memory[["submodels"]])) {
      for (position in seq_along(in_memory[["submodels"]][[submodel]])) {
        expect_equal(draws_of(back, submodel, position),
                     draws_of(in_memory, submodel, position),
                     info = paste(type, submodel, position))
        expect_equal(back[["submodels"]][[submodel]][[position]][["posterior"]][["loglik"]],
                     in_memory[["submodels"]][[submodel]][[position]][["posterior"]][["loglik"]],
                     info = paste(type, submodel, position))
      }
    }
  }
})

test_that("the draws do not depend on the number of workers", {
  skip_on_cran()

  object <- export_model("gvar")
  one <- written_export(object, "cores-one")
  two <- written_export(object, "cores-two")

  estimate <- function(export, cores) {
    export <- add_priors(export, coef = list(v_i = 1),
                         sigma = list(df = 3, scale = 1e-8), cores = cores)
    export <- add_initial_values(export, cores = cores)
    export <- add_seed(export, 99, cores = cores)
    add_posterior_coefficients(export, cores = cores)
  }
  one <- estimate(one, 1)
  two <- estimate(two, 2)

  first <- read_gvar_from_folder(one[["folder"]])
  second <- read_gvar_from_folder(two[["folder"]])
  for (submodel in names(first[["submodels"]])) {
    expect_equal(draws_of(second, submodel, 1), draws_of(first, submodel, 1),
                 info = submodel)
  }
})

test_that("the manifest follows what the models become", {
  object <- export_model("gvar")
  export <- written_export(object, "manifest")
  export <- add_priors(export, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))
  export <- add_posterior_coefficients(add_initial_values(export))

  before <- export[["manifest"]]

  # Thinning takes draws out of the files. The manifest describes what the
  # sampler was asked for rather than what is left of its draws, so it stays as
  # it was -- but it is written back, and reading the folder agrees with it.
  export <- thin(export, thin = 2)
  expect_equal(export[["manifest"]], before)
  expect_equal(read_export(export[["folder"]])[["manifest"]], before)

  drawn <- read_gvar_from_folder(export[["folder"]])
  expect_equal(nrow(draws_of(drawn, names(drawn[["submodels"]])[1], 1)), 10L)

  # A change of the specification itself does reach the manifest.
  changed <- map_submodels(export, function(model) {
    model[["model"]][["burnin"]] <- 77L
    model
  })
  expect_true(all(changed[["manifest"]][["burnin"]] == 77L))
  expect_equal(read_export(export[["folder"]])[["manifest"]], changed[["manifest"]])
})

test_that("map_submodels passes the position and refuses what is not a model", {
  object <- export_model("gvar")
  export <- written_export(object, "map")

  seen <- map_submodels(export, function(model, index) {
    model[["model"]][["iterations"]] <- as.integer(100 + index)
    model
  })
  expect_identical(seen[["manifest"]][["iterations"]],
                   as.integer(100 + seq_len(nrow(seen[["manifest"]]))))

  expect_error(map_submodels(export, function(model) "not a model"),
               "rather than a model")
  expect_error(map_submodels(object, function(model) model), "'gvarexport'")
})

test_that("a model on disk is left alone when a step fails", {
  object <- export_model("gvar")
  export <- written_export(object, "failing")
  file <- file.path(export[["folder"]], export[["manifest"]][1, "file"])
  before <- bayests::read_model_from_hdf5(filename = file)

  expect_error(map_submodels(export, function(model) stop("no")), "no")
  after <- bayests::read_model_from_hdf5(filename = file)
  expect_equal(after[["data"]][["train"]][["y"]], before[["data"]][["train"]][["y"]])
  expect_false(file.exists(paste0(file, ".writing")))
})

# --- Solving from the export ------------------------------------------------------
#
# A global model stacked from pieces of the chain has to be the model stacked
# from the whole of it, draw for draw.

solved_pair <- function(type, chunk) {
  object <- export_model(type)
  export <- written_export(object, paste0("solve-", type, "-", chunk))

  in_memory <- export_priors(object, type)
  in_memory <- add_seed(add_initial_values(in_memory), 4321)
  in_memory <- add_posterior_coefficients(in_memory)

  on_disk <- add_seed(add_initial_values(export_priors(export, type)), 4321)
  on_disk <- add_posterior_coefficients(on_disk)

  positions <- if (type == "gvec") {
    manifest <- on_disk[["manifest"]]
    do.call(rbind, lapply(unique(manifest[["submodel"]]), function(submodel) {
      data.frame(submodel = submodel, position = 1L)
    }))
  } else {
    NULL
  }
  if (type == "gvec") {
    # One model per sub-model, as stacking requires.
    in_memory[["submodels"]] <- lapply(in_memory[["submodels"]], function(models) {
      models[1]
    })
  }

  list(memory = suppressMessages(submodels_to_gvar(in_memory)),
       disk = suppressMessages(submodels_to_gvar(on_disk, submodels = positions,
                                                 chunk = chunk)))
}

test_that("a global model solved from pieces is the one solved whole", {
  for (type in c("gvar", "gvec")) {
    for (chunk in c(7, 1000)) {
      pair <- solved_pair(type, chunk)
      expect_equal(pair[["disk"]][["posterior"]][["a"]][["coeffs"]],
                   pair[["memory"]][["posterior"]][["a"]][["coeffs"]],
                   info = paste(type, chunk))
      expect_equal(pair[["disk"]][["posterior"]][["u_sigma_inv"]][["coeffs"]],
                   pair[["memory"]][["posterior"]][["u_sigma_inv"]][["coeffs"]],
                   info = paste(type, chunk))
      expect_equal(pair[["disk"]][["g"]], pair[["memory"]][["g"]],
                   info = paste(type, chunk))
      expect_equal(pair[["disk"]][["data"]], pair[["memory"]][["data"]],
                   info = paste(type, chunk))
      expect_identical(pair[["disk"]][["model"]], pair[["memory"]][["model"]],
                       info = paste(type, chunk))
      # The seed included, which both paths take from the sub-models' seeds.
      expect_false(is.null(pair[["disk"]][["model"]][["seed"]]), info = paste(type, chunk))
      expect_s3_class(pair[["disk"]], "bvarmodel")
    }
  }
})

test_that("an export with several models per sub-model says which to name", {
  object <- export_model("gvec")
  export <- written_export(object, "solve-ambiguous")
  export <- add_posterior_coefficients(add_initial_values(export_priors(export, "gvec")))
  expect_error(submodels_to_gvar(export), "more than one model")
})

test_that("a step over one sub-model at a time is the step over all of them", {
  skip_on_cran()

  object <- export_model("gvar")
  whole <- written_export(object, "subset-whole")
  piecewise <- written_export(object, "subset-piecewise")

  estimate <- function(export, submodels) {
    export <- add_priors(export, coef = list(v_i = 1),
                         sigma = list(df = 3, scale = 1e-8),
                         submodels = submodels)
    export <- add_initial_values(export, submodels = submodels)
    export <- add_seed(export, 99, submodels = submodels)
    add_posterior_coefficients(export, submodels = submodels)
  }

  whole <- estimate(whole, NULL)
  for (submodel in unique(piecewise[["manifest"]][, "submodel"])) {
    piecewise <- estimate(piecewise, submodel)
  }

  first <- read_gvar_from_folder(whole[["folder"]])
  second <- read_gvar_from_folder(piecewise[["folder"]])
  for (submodel in names(first[["submodels"]])) {
    expect_equal(draws_of(second, submodel, 1), draws_of(first, submodel, 1),
                 info = submodel)
  }
  expect_equal(read_export(piecewise[["folder"]])[["manifest"]],
               read_export(whole[["folder"]])[["manifest"]])
})

test_that("a step leaves the sub-models it was not asked for alone", {

  object <- export_model("gvar")
  export <- written_export(object, "subset-untouched")
  submodels <- unique(export[["manifest"]][, "submodel"])

  before <- vapply(export[["manifest"]][, "file"], function(file) {
    file.info(file.path(export[["folder"]], file))[["mtime"]]
  }, numeric(1))

  Sys.sleep(1)
  export <- add_priors(export, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 1e-8),
                       submodels = submodels[1])

  after <- vapply(export[["manifest"]][, "file"], function(file) {
    file.info(file.path(export[["folder"]], file))[["mtime"]]
  }, numeric(1))

  touched <- export[["manifest"]][, "submodel"] == submodels[1]
  expect_true(all(after[touched] > before[touched]))
  expect_equal(after[!touched], before[!touched])

  # The models of the other sub-models have no priors, so they cannot be
  # estimated, while the one that was asked for can.
  models <- read_gvar_from_folder(export[["folder"]])
  expect_false(is.null(models[["submodels"]][[submodels[1]]][[1]][["priors"]]))
  expect_null(models[["submodels"]][[submodels[2]]][[1]][["priors"]])

  expect_error(map_submodels(export, add_initial_values, submodels = "XX"),
               "does not contain the sub-model")
})


test_that("an export names the directory of each of its sub-models", {

  object <- export_model("gvec")
  export <- written_export(object, "directories")
  submodels <- unique(export[["manifest"]][, "submodel"])

  directories <- submodel_directories(export)
  expect_identical(names(directories), submodels)
  expect_true(all(dir.exists(directories)))
  for (submodel in submodels) {
    files <- list.files(directories[[submodel]], pattern = "[.]h5$")
    expect_equal(length(files),
                 sum(export[["manifest"]][, "submodel"] == submodel),
                 info = submodel)
  }

  one <- submodel_directories(export, submodels = submodels[1])
  expect_identical(one, directories[submodels[1]])

  expect_error(submodel_directories(export, submodels = "XX"),
               "does not contain the sub-model")
  expect_error(submodel_directories(unclass(export)), "must be of class")
})
