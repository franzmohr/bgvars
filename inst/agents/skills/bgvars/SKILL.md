---
name: bgvars
description: Write correct bgvars code — Bayesian global VAR (GVAR) and global VEC (GVEC) models in R with create_gvarmodel, create_gvecmodel, add_weight_matrices, add_submodels, create_varxsubmodel, create_vecxsubmodel, align_model_obs, add_priors, add_posterior_coefficients, scale_error_correction, write_to_hdf5, open_export, map_submodels, choose_best_model_from_hdf5, read_gvar_from_folder, submodels_to_gvar and persistence_profiles. Use whenever R code loads bgvars, handles a 'gvarmodel', 'gvecmodel', 'submodeldata', 'varxsubmodel', 'vecxsubmodel' or 'gvarexport' object, uses the gvar2019, gvar2023 or dees2007 data, or estimates country models with weakly exogenous foreign variables and trade weights. Also use when a global model is estimated through a folder of HDF5 files, in batches, or by the BayesTS command line, and when the solved global model is analysed with generalised impulse responses.
---

# Writing correct bgvars code

bgvars estimates Bayesian **global** VAR and VEC models: one VARX\* or VECX\*
model per country (a *sub-model*), each taking weighted averages of the other
countries' series as weakly exogenous regressors, stacked afterwards into one
VAR for all countries. It is built on bvartools. The sub-models **are**
bvartools models (a `varxsubmodel` is a `bvarmodel`, a `vecxsubmodel` a
`bvecmodel`), the verbs are bvartools generics, and the solved global model is a
plain `bvarmodel`. The Gibbs samplers are the C++ core of BayesTS.

```r
library(bgvars)
set.seed(42)

data("gvar2023")
submodel_data <- select_list_elements(gvar2023[["submodel_data"]], c("DE", "JP", "US"))
global_data <- gvar2023[["global_data"]]

object <- create_gvarmodel(submodel_data = submodel_data, global_data = global_data)
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
object <- add_submodels(object,
                        endogen = c("y", "Dp", "r"), p_endogen = 1,
                        exogen = c("y", "Dp", "r"), p_exogen = 2,
                        global = "poil", s = 1,
                        iterations = 200, burnin = 100)
object <- align_model_obs(object)
object <- add_priors(object,
                     coef = list(v_i = 0),
                     sigma = list(df = 3, scale = 1e-4))
object <- add_initial_values(object)
object <- add_posterior_coefficients(object)

gvar <- submodels_to_gvar(object)
gir <- irf(gvar, impulse = "US_r", response = "DE_y", n_ahead = 8,
           type = "gir", shock = "sd")

stopifnot(inherits(gvar, "bvarmodel"), gvar[["model"]][["k"]] == 9)
```

`iterations = 200` keeps the example quick. The defaults are 10000 draws after
2000 burn-in per sub-model; `thin = t` runs the chain `t` times as long and keeps
the same number of draws.

## The rules that prevent wrong results

**1. Assign every step back.** Each function returns the object with something
added: `object <- add_priors(object, ...)`, never `add_priors(object, ...)` alone.

**2. The order is fixed.** `create_gvarmodel()` or `create_gvecmodel()`, then
`add_weight_matrices()`, then the sub-models (`add_submodels()` for every
country at once, or `create_varxsubmodel()` / `create_vecxsubmodel()` per
country), then optionally `align_model_obs()`, then `add_priors()`,
`add_initial_values()`, `add_posterior_coefficients()`, and
`add_posterior_loglik()` when specifications are compared. Weights come before
sub-models because a sub-model's foreign variables are built from them.

**3. `add_priors()` has no defaults.** `coef` and `sigma` are always required,
`coint` for a GVEC. Prior variances are **precisions**: `v_i = 0` is
uninformative. Which `sigma` elements are needed depends on `error`, exactly as
in bvartools (`df` and `scale` for the Wishart default, `shape` and `rate` for
gamma, six elements for stochastic volatility). Time varying sub-models also need
a prior for their state variances, in one of two parameterisations: `shape` and
`rate` in `coef` put a gamma prior on the precisions of the state equation, and
`omega_v` puts a normal prior on the signed standard deviations instead, the
non-centred form of Fruehwirth-Schnatter and Wagner (2010). The second makes a
constant coefficient an interior point of the prior rather than the boundary of a
variance, so the posterior can say that a coefficient does not move: it carries
the draws of omega and, in `omega_log_zero`, the Savage-Dickey ordinate that
tests omega = 0 per state. `omega_v` is also the prior mean of the state
variance, which is `rate / (shape - 1)` for the other, so the two are matched by
setting them equal. They cannot be combined, a VECX sub-model takes
`omega_v_alpha` for its loadings as it takes `rate_alpha`, and an unrecognised
element stops with an error. (Some help pages still speak of default values; the
code has none.)

**4. Lag orders of weakly exogenous and global variables count from lag
zero.** `p_exogen = 1` means the contemporaneous foreign variables alone, `2`
adds their first lag, `0` leaves them out; `s` counts the global variables the
same way. For a VECX sub-model `p_endogen`, `p_exogen` and `s` are lag orders
**of the model in levels**, so a VARX\*(p, q) in the notation of Dees et al.
(2007) is `p_endogen = p, p_exogen = q + 1`, and `p_endogen = 1` has no lagged
differences at all.

**5. Vectors make candidates.** A vector for `p_endogen`, `p_exogen`, `s`, `r`,
`delta_beta` or `delta_sigma` gives every sub-model one candidate per
combination, stored as a `modellist` in `object$submodels[[country]]`. Every step
maps over all of them. **`submodels_to_gvar()` needs exactly one model per
sub-model** and stops otherwise; reduce the candidates first (rule 8).

**6. The solved model is the reduced form.** `submodels_to_gvar()` stacks the
sub-models into `G y_t = A_1 y_{t-1} + ... + u_t` and returns
`y_t = G^{-1} A_1 y_{t-1} + ... + G^{-1} u_t` as a `bvarmodel`, draw by draw, with
the draws of `G` in element `g`. Its variables are named `<country>_<variable>`
(`US_r`, `DE_y`). A GVEC is solved in levels, each sub-model rewritten by
`vec_to_var()` first. Everything bvartools does with a `bvarmodel` —
`irf(type = "gir")`, `fevd()`, `add_forecast_input()` with the future path of
the global variables in `exogen` — works on it.

A GVEC keeps one thing more than the level form holds. `vec_to_var()` drops the
cointegrating vectors, so `submodels_to_gvar()` takes them first and maps each
through the weight matrix of its sub-model: the solved model carries them in
`posterior$beta`, counts them in `model$rank`, and is classed `gvecsolved`
before `bvarmodel`. `persistence_profiles(gvar)` then follows each relation
after a system-wide shock — one on impact, falling to zero where the relation is
cointegrating, and the check on whether the rank was real. The same function on
a single sub-model warns that it can follow only the domestic part, since a
sub-model's foreign variables are not endogenous to it. A sub-model already in
level form contributes no vectors, so a model stacked from both forms carries
fewer relations than its ranks suggest and `model$rank` counts the ones it has.

**7. Time variation is solved one period at a time.** A sub-model with
`tvp = TRUE` or stochastic volatility has one set of coefficients per period;
`submodels_to_gvar(object, period = ...)` solves the global model at that date
(default: the last period of the sample). **`period` is a row number of the
global data, not a time** like `2019.75`. Call it once per period of interest.

**8. Compare specifications per sub-model.** `add_posterior_loglik()`, then
`selection_criteria()` on one sub-model's `modellist`, or
`choose_best_model_from_hdf5(folder, criterion = "BIC")` for all of them from an
export. That returns a data frame with columns `submodel` and `position`, which
`read_gvar_from_folder(folder, submodels = best)` and
`submodels_to_gvar(export, submodels = best)` take directly. Use a common
sample (`align_model_obs()`) before comparing lag orders.

**9. Scaling the error correction term is a modelling step.** For a GVEC with a
prior on the cointegration space, `scale_error_correction(object)` before the
simulation puts the series of the error correction term on comparable scales,
and `rescale_error_correction(object)` afterwards puts series and draws back.
**Rescale before `vec_to_var()` or `submodels_to_gvar()`.** The estimation never
scales on its own.

**10. Seeds.** `set.seed()` reaches the samplers. `add_seed(object, seed)` gives
the models `seed`, `seed + 1`, … in order, so the draws do not depend on how many
`cores` a step runs on or whether it ran in one go or several.

## Large models: work through a folder

A global model is large in its draws, and a time varying sub-model keeps a
coefficient path per draw — hundreds of megabytes each. The sub-models are
independent given the weights, so bgvars can hold **one model at a time** by
working on an export rather than on an object. This is the standard workflow for
anything beyond a small constant-coefficient model, not a fallback:

```r
folder <- file.path(tempdir(), "skill-gvar")
dir.create(folder, showWarnings = FALSE)
write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)

export <- open_export(folder)
export <- add_seed(export, 20260919)
export <- add_posterior_coefficients(export)
export <- add_posterior_loglik(export)

best <- choose_best_model_from_hdf5(folder, criterion = "BIC")
gvar <- submodels_to_gvar(export, submodels = best)

stopifnot(inherits(export, "gvarexport"), setequal(best$submodel, c("DE", "JP", "US")))
```

The export is `model.h5` (global data, weights, a manifest) and one file per
model under `submodels/<country>/<nnn>.h5`. Every estimation verb has a method
for the `gvarexport` handle, with `cores` and a `submodels` argument that
restricts a step to some countries (which is how an interrupted run resumes).
`map_submodels(export, f)` applies any other function model by model.
`write_to_hdf5(..., append = TRUE)` adds batches of sub-models to a folder, so a
model too large to set up at once is built batch by batch. The files are BayesTS
model files, so `bayests posterior <folder>/submodels` estimates them outside R.
Details and the batch loop: `references/export.md`.

## Where to read next

| File | Contents |
| --- | --- |
| `references/recipes.md` | Complete examples: a GVAR with candidate lag orders chosen in memory, the Dees et al. (2007) GVEC with oil endogenous in the US model and scaled error correction terms, one sub-model on its own, a TVP-SV GVAR through a folder |
| `references/specification.md` | The data (`submodeldata`), weights, what `add_submodels()` builds, naming of foreign variables, lag counting, deterministic terms, priors per model type, what is refused |
| `references/export.md` | The folder workflow: batches, `open_export()`, `map_submodels()`, BayesTS, reading back, choosing and solving from files |
| `references/solving.md` | `submodels_to_gvar()`, the layout of the solved model, `period`, generalised impulse responses and what they shock, stability |

The help pages are the reference; this skill summarises them. Read the method,
not the generic: `?add_priors.varxsubmodel`, `?add_priors.vecxsubmodel`,
`?create_vecxsubmodel`, `?submodels_to_gvar`, `?write_to_hdf5.gvarmodel`,
`?open_export`. The vignettes (`browseVignettes("bgvars")`) cover a GVAR, a GVEC
replicating Dees et al. (2007), and TVP-SV versions of both. The copy of this
skill matching the installed version is at
`system.file("agents", package = "bgvars")`.

## Related skills

- **bvartools** — everything about a single model: the prior lists in full, the
  layout of `posterior`, `irf()`, `fevd()`, forecasts, `selection_criteria()`.
  A sub-model is a bvartools model, so that skill applies to it unchanged.
- **bayests** — the HDF5 model file format and the `bayests` command line, for
  when the sub-model files of an export are estimated or inspected outside R.
