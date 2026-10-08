# Estimating through a folder

A global model is large in its posterior draws, not in its data. A sub-model
with time varying coefficients keeps a coefficient path per draw — some 0.6 GB
for six variables at 3000 draws — so thirty of them need tens of gigabytes in
one session. Given the weight matrices the sub-models are independent, so
bgvars can hold one model at a time by working on files. For anything beyond a
small model with constant coefficients, this is the workflow, and the in-memory
object is the special case.

Every `r` block here runs in the package's test suite.

## The export

`write_to_hdf5(object, folder, overwrite = FALSE, append = FALSE, mc.cores)`
writes

| Path | Holds |
| --- | --- |
| `model.h5` | the global data, one weight matrix per sub-model, and the **manifest**: one row per model file with its sub-model, algorithm, lag orders, rank, chain length |
| `submodels/<country>/<nnn>.h5` | one estimable model each — data, priors, starting values, and later its draws — numbered in the order `add_submodels()` produced the candidates |

Each sub-model file is a complete BayesTS model file (see the `bayests`
skill). `overwrite = FALSE` stops rather than touch an existing export.

```r
library(bgvars)
set.seed(1)
data("gvar2023")
submodel_data <- select_list_elements(gvar2023[["submodel_data"]], c("DE", "JP", "US"))

object <- create_gvarmodel(submodel_data, global_data = gvar2023[["global_data"]])
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
object <- add_submodels(object,
                        endogen = c("y", "Dp", "r"), p_endogen = 1:2,
                        exogen = c("y", "Dp", "r"), p_exogen = 1,
                        iterations = 100, burnin = 50)
object <- align_model_obs(object)

folder <- file.path(tempdir(), "skill-export")
dir.create(folder, showWarnings = FALSE)
write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)

stopifnot(file.exists(file.path(folder, "model.h5")),
          length(list.files(file.path(folder, "submodels"), recursive = TRUE)) == 6)
```

## Estimating the export in R

`read_export(folder)` returns a `gvarexport` handle: the folder, the manifest and
the class the export was written from — not the models. Every estimation verb has
a method for it that reads a model, changes it, writes it back (through a
temporary file, so an interruption never leaves half a model) and drops it:
`add_priors()`, `add_initial_values()`, `add_seed()`,
`add_posterior_coefficients()`, `add_posterior_loglik()`, `thin()`,
`scale_error_correction()` and `rescale_error_correction()`.

```r
export <- read_export(folder)
export <- add_priors(export, coef = list(v_i = 0), sigma = list(df = 3, scale = 1e-4))
export <- add_initial_values(export)
export <- add_seed(export, 20260919)
export <- add_posterior_coefficients(export, cores = 1)
export <- add_posterior_loglik(export)

stopifnot(inherits(export, "gvarexport"), nrow(export$manifest) == 6)
```

- `cores` spreads the models over a socket cluster (works on Windows). Each
  model has its own seed, so the draws do not depend on `cores`.
- `submodels = c("DE", "US")` restricts a step to those countries and leaves the
  rest untouched. Seeds are numbered by the position in the **whole** manifest,
  so estimating country by country draws exactly what one run draws. That is
  how an interrupted estimation is resumed.
- `map_submodels(export, f, ..., submodels, cores)` is the general form: `f`
  takes a single `bvarmodel` or `bvecmodel` and returns one. If `f` has an
  argument `index`, it gets the model's position in the manifest. With several
  cores `f` is sent to workers that load the *installed* packages — pass
  anything it needs as an argument or through `export =`.

```r
export <- map_submodels(export, function(model, index) {
  # index is the model's row in the manifest
  stopifnot(index %in% 3:4, !is.null(model[["posterior"]]))
  model
}, submodels = "JP")
```

## Estimating the export outside R

The sub-model files are what BayesTS reads, so the whole folder is estimated
by the command line, which writes the draws back into the files it read:

```text
bayests posterior <folder>/submodels
```

From R, `submodel_directories(export)` gives one directory per country, and
`bayests::bayests_files()` runs several BayesTS processes over them — the
parallelism is in processes, not in R workers, which on Windows cost gigabytes
each while they only wait:

```text
run <- bayests::bayests_files(executable = "D:/BayesTS/bin/bayests.exe")
run(submodel_directories(export), jobs = 6)
```

`bayests::bayests_posterior()` wraps the same executable as a
`posterior_function` for `add_posterior_coefficients()`. A sub-model with
`algorithm = "discount"` can only be estimated this way.

## Building the folder in batches

When even the unestimated sub-models are more than a session should hold, or a
country is added to an estimated model, write batches with `append = TRUE`. The
first write creates the export; later ones add their sub-models' files (or
replace files already there) and merge their manifest rows, leaving global data
and weights as the first write left them. Estimating inside the loop, before
each write, is equally valid.

```r
submodels <- names(object[["submodels"]])
batches <- split(submodels, rep(1:2, length.out = length(submodels)))
batch_folder <- file.path(tempdir(), "skill-batches")
dir.create(batch_folder, showWarnings = FALSE)

for (i in seq_along(batches)) {
  batch <- object
  batch[["submodels"]] <- object[["submodels"]][batches[[i]]]
  write_to_hdf5(batch, folder = batch_folder, overwrite = i == 1,
                append = i > 1, mc.cores = 1)
  rm(batch)
}

stopifnot(nrow(read_export(batch_folder)$manifest) == 6)
```

## Reading back, choosing, solving

Nothing downstream needs the whole model in memory:

| Function | Reads |
| --- | --- |
| `choose_best_model_from_hdf5(folder, criterion)` | only the log-likelihood draws for `"LL"`, `"WAIC"`, `"LOOIC"`; each file in full, one at a time, for `"AIC"`, `"BIC"`, `"HQ"`. Returns `data.frame(submodel, position)`. A model without draws is skipped with a warning; a sub-model with none gets `NA` |
| `read_submodel_from_folder(folder, submodel)` | the candidates of one country as a `modellist`, for `selection_criteria()` |
| `read_gvar_from_folder(folder, submodels = best, draws = NULL)` | the global model with one model per sub-model, draws and scaling factors included, found through the manifest (so a moved or renamed folder still works) |
| `submodels_to_gvar(export, submodels = best, chunk = 100)` | the chosen models in pieces of `chunk` draws, so it never holds every chain at once |

```r
best <- choose_best_model_from_hdf5(folder, criterion = "BIC")
candidates <- read_submodel_from_folder(folder, submodel = "US")
sc <- selection_criteria(candidates)

gvar <- submodels_to_gvar(export, submodels = best)

stopifnot(names(best) == c("submodel", "position"),
          inherits(candidates, "modellist"), length(candidates) == 2,
          inherits(gvar, "bvarmodel"))
```

The solved model is held whole. For a large one, write it with
`write_to_hdf5(gvar, filename = file)` and analyse it through
`bayests::read_model_file(file)`, whose `irf()` and `fevd()` read the draws in
pieces too.

## Errors worth recognising

- *"may only contain one model per sub-model"* — pass `submodels = best` or
  choose candidates first.
- *"carry a closed-form posterior rather than draws"* — a discounted sub-model
  was chosen for solving; re-estimate that country with a sampler.
- A step with `cores > 1` failing on a function that works in the session — the
  workers load the installed package, not `devtools::load_all()`; install first.
