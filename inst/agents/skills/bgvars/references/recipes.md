# Recipes

Complete, runnable analyses. Every `r` block runs in the package's test suite, so
the shapes asserted with `stopifnot()` are what the installed version produces.
The chains are far too short for inference; use thousands of draws in practice.

## A GVAR with candidate lag orders, chosen in memory

Each country gets two candidates. Comparing them needs the log-likelihood and a
common sample; the chosen candidate then replaces the `modellist` of its country.

```r
library(bgvars)
set.seed(11)
data("gvar2023")
submodel_data <- select_list_elements(gvar2023[["submodel_data"]], c("DE", "JP", "US"))

object <- create_gvarmodel(submodel_data, global_data = gvar2023[["global_data"]])
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
object <- add_submodels(object,
                        endogen = c("y", "Dp", "r"), p_endogen = 1:2,
                        exogen = c("y", "Dp", "r"), p_exogen = 1,
                        global = "poil", s = 1,
                        iterations = 200, burnin = 100)
object <- align_model_obs(object)
object <- add_priors(object, coef = list(v_i = 0), sigma = list(df = 3, scale = 1e-4))
object <- add_initial_values(object)
object <- add_seed(object, 20260924)
object <- add_posterior_coefficients(object)
object <- add_posterior_loglik(object)

for (country in names(object[["submodels"]])) {
  sc <- selection_criteria(object[["submodels"]][[country]])
  best <- choose_best_model(sc, criterion = "BIC")[1]   # ties return several
  object[["submodels"]][[country]] <- object[["submodels"]][[country]][best]
}

gvar <- submodels_to_gvar(object)
stopifnot(all(lengths(object[["submodels"]]) == 1), inherits(gvar, "bvarmodel"))
```

`object$submodels[[country]][best]` (single brackets) keeps the `modellist`
class; `[[best]]` would drop it. For an export, `choose_best_model_from_hdf5()`
does this loop without reading the coefficient draws (see `export.md`).

## The Dees et al. (2007) GVEC: oil endogenous in the US model

The oil price is added to the US series, is endogenous there, and enters every
other country as a weakly exogenous variable. Sub-models are VECX\*(p, q) with
an unrestricted constant and a trend restricted to the cointegration space,
estimated with flat priors on scaled error correction terms.

```r
library(bgvars)
set.seed(2007)
data("dees2007")
submodel_data <- select_list_elements(dees2007[["submodel_data"]], c("EA", "JP", "US"))
global_data <- dees2007[["global_data"]]

us <- submodel_data[["US"]][["endogen"]]
vars_us <- c(dimnames(us)[[2]], "poil")
submodel_data[["US"]][["endogen"]] <- cbind(us, global_data)
dimnames(submodel_data[["US"]][["endogen"]])[[2]] <- vars_us

object <- create_gvecmodel(submodel_data = submodel_data)
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 1999:2001)

spec <- data.frame(submodel = c("EA", "JP", "US"), p = c(2, 1, 2), q = c(2, 1, 2), r = c(2, 1, 2))
for (i in seq_len(nrow(spec))) {
  country <- spec$submodel[i]
  us_model <- country == "US"
  object[["submodels"]][[country]] <- create_vecxsubmodel(object,
    submodel = country,
    endogen = if (us_model) c("y", "Dp", "rs", "rl", "poil") else c("y", "Dp", "rs", "rl", "ep"),
    p_endogen = spec$p[i],
    exogen = if (us_model) c("y", "Dp") else c("y", "Dp", "rs", "poil"),
    p_exogen = spec$q[i] + 1,
    r = spec$r[i],
    const = "unrestricted", trend = "restricted",
    iterations = 200, burnin = 100)
}
object <- align_model_obs(object)

object <- add_priors(object,
                     coef = list(v_i = 0),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 3, scale = 1e-4))
object <- add_initial_values(object)
object <- scale_error_correction(object)
object <- add_posterior_coefficients(object)
object <- rescale_error_correction(object)

gvec <- submodels_to_gvar(object)
stopifnot(inherits(gvec, "bvarmodel"), "US_poil" %in% gvec[["model"]][["endogen"]],
          !"EA_poil" %in% gvec[["model"]][["endogen"]])
```

`add_initial_values()` starts each chain at the reduced-rank ML estimate of its
sub-model, which is also what the paper reports. With flat priors the posterior
averages over a diffuse posterior of the cointegration space, so many drawn
global models may be unstable; see `solving.md`.

## One sub-model on its own

For choosing a country's lag order or rank before building the global model,
work on its `modellist` directly. `r = NULL` gives candidates for every rank.

```r
library(bgvars)
set.seed(3)
data("dees2007")
submodel_data <- select_list_elements(dees2007[["submodel_data"]], c("EA", "JP", "US"))

object <- create_gvecmodel(submodel_data = submodel_data)
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 1999:2001)

ea <- create_vecxsubmodel(object, submodel = "EA",
                          endogen = c("y", "Dp", "rs"), p_endogen = 2,
                          exogen = c("y", "Dp", "rs"), p_exogen = 2,
                          r = NULL, const = "unrestricted", trend = "restricted",
                          iterations = 200, burnin = 100)
ea <- add_priors(ea, coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
                 sigma = list(df = 3, scale = 1e-4))
ea <- add_initial_values(ea)
ea <- add_posterior_coefficients(ea)
ea <- add_posterior_loglik(ea)
sc <- selection_criteria(ea)

stopifnot(inherits(ea, "modellist"), length(ea) == 4)   # ranks 0 to 3
```

A `vecxsubmodel` is a `bvecmodel`, so `summary()`, `plot()` and
`vec_to_var()` work on a single element such as `ea[[2]]`.

## A time varying GVAR with stochastic volatility, through a folder

The configuration of the TVP-SV-GVAR vignette, sized down. Such models are
estimated through an export (see `export.md`); `period` picks the date at which
the global model is solved, as a **row number of the global data**, not a time.

```r
library(bgvars)
data("gvar2023")
submodel_data <- select_list_elements(gvar2023[["submodel_data"]], c("DE", "US"))
object <- create_gvarmodel(submodel_data, global_data = gvar2023[["global_data"]])
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
object <- add_submodels(object,
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 1,
                        tvp = TRUE, error = "sv", varsel = "bvs",
                        iterations = 30, burnin = 20)

folder <- file.path(tempdir(), "skill-tvp")
dir.create(folder, showWarnings = FALSE)
write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)

export <- open_export(folder)
export <- add_priors(export,
                     coef = list(v_i = 1, v_i_det = 1 / 10, shape = 3, rate = 1e-10, rate_det = 1e-8),
                     sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                  state_variance = 0.05, offset = 1e-8),
                     varsel = list(inprior = 0.5, exclude_det = TRUE))
export <- add_initial_values(export)
export <- add_seed(export, 1)
export <- add_posterior_coefficients(export)

global_time <- stats::time(object[["global"]][["endogen"]])
period_2019q4 <- which(abs(global_time - 2019.75) < 1e-6)
gvar_2019 <- submodels_to_gvar(export, period = period_2019q4)
stopifnot(inherits(gvar_2019, "bvarmodel"), gvar_2019[["model"]][["k"]] == 4)
```
