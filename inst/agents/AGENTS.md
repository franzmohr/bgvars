# AGENTS.md — using bgvars

Guidance for coding agents writing R code with
[bgvars](https://github.com/franzmohr/bgvars), which estimates Bayesian global
VAR (GVAR) and global VEC (GVEC) models on top of
[bvartools](https://github.com/franzmohr/bvartools). The full skill is
`skills/bgvars/SKILL.md`, with tested references in `skills/bgvars/references/`:
`recipes.md` for complete examples, `specification.md` for data, weights,
sub-models and priors, `export.md` for estimating through a folder of HDF5
files, and `solving.md` for the solved global model and its analysis. This file
is the part worth keeping in context at all times.

1. **Assign every step back**: `object <- add_priors(object, ...)`.
2. **The order**: `create_gvarmodel()` or `create_gvecmodel()`, then
   `add_weight_matrices()`, then `add_submodels()` (or `create_varxsubmodel()` /
   `create_vecxsubmodel()` per country, assigned into `object$submodels`), then
   `align_model_obs()`, `add_priors()`, `add_initial_values()`,
   `add_posterior_coefficients()`, and `add_posterior_loglik()` for comparisons.
   `submodels_to_gvar()` then solves the global model.
3. **`add_priors()` has no defaults**: `coef` and `sigma` always, `coint` for a
   GVEC. Prior variances are precisions; `v_i = 0` is flat.
4. **Weakly exogenous and global lags count from lag zero**: `p_exogen = 1` is
   the contemporaneous term alone. VECX lag orders are those of the model in
   levels: VARX\*(p, q) is `p_endogen = p, p_exogen = q + 1`.
5. **Vectors make candidates** per country; `submodels_to_gvar()` needs exactly
   one model per sub-model. Choose with `selection_criteria()` or
   `choose_best_model_from_hdf5()`, whose `data.frame(submodel, position)` goes
   to `read_gvar_from_folder()` and `submodels_to_gvar()` as `submodels`.
6. **The solved model is a reduced-form `bvarmodel`** with variables named
   `<country>_<variable>` and the draws of `G` in `$g`. Use bvartools on it:
   `irf(type = "gir")`, `fevd()`, forecasts with the global series' future path
   in `exogen`. Time varying sub-models are solved at one `period` at a time.
7. **GVEC: scale, estimate, rescale**: `scale_error_correction()` before the
   simulation, `rescale_error_correction()` before `vec_to_var()` or
   `submodels_to_gvar()`.
8. **Large models go through a folder**: `write_to_hdf5(object, folder)` (in
   batches with `append = TRUE`), `open_export(folder)`, then the same verbs on
   the handle, with `cores` and `submodels`. The files are BayesTS model files:
   `bayests posterior <folder>/submodels` estimates them outside R.
9. **Seeds**: `add_seed(object, seed)` numbers the models, so draws do not
   depend on `cores` or on how a run was split.
10. **Refusals are statistical**: no SSVS for VECX sub-models, no structural
    sub-model with an unrestricted covariance, no global model from discounted
    (closed-form) sub-models.

The installed package carries these files at
`system.file("agents", package = "bgvars")`, matching its version. For the
internals of a single sub-model — prior lists, posterior layout, `irf()` — the
bvartools guide applies unchanged.
