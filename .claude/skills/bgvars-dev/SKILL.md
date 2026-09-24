---
name: bgvars-dev
description: How to change the bgvars R package itself -- which layer a change belongs in (BayesTS, bvartools or bgvars), how a verb is laid out across the global model, the export handle and the sub-model classes, the roxygen, test and vignette conventions, the commands that build, test and check it on this machine, the JSS manuscript, and the branch-and-merge workflow. Use for any development work inside the bgvars repository -- adding or changing an exported function or method, fixing a bug, writing tests, touching the HDF5 export, regenerating documentation or vignettes, preparing a release, or reviewing a diff -- even when the request does not say "develop". For writing analysis code that merely *uses* bgvars, use the bgvars skill instead.
---

# Developing bgvars

bgvars estimates Bayesian global VAR and VEC models. It is the top of a stack
of three of the user's own projects, and most design questions are answered by
which layer a piece of work belongs to:

| Layer | Owns |
| --- | --- |
| **BayesTS** (`D:/workspace-cpp/BayesTS`, C++) | the samplers and everything computationally heavy; the HDF5 model file format and the `bayests` command line |
| **bvartools** (`D:/workspace-r/bvartools`) | a single VAR/VEC model in R: the model object, priors, the R bindings to the vendored samplers, `irf()`, `fevd()`, forecasts, `selection_criteria()`, `write_to_hdf5()` for one model, `open_model()`, `bayests_files()` |
| **bgvars** (this repository) | what only a *global* model has: country data and trade weights, foreign variables, sub-models with two blocks of unmodelled variables, the export of many models to a folder, choosing per sub-model, and stacking the sub-models into the solved global model |

This skill is about changing the package. The user-facing skill that teaches
how to *use* it lives in `inst/agents/skills/bgvars/` and ships with the
package; it is one of the things a change may have to update (see *Finishing a
change*). The `bayests` and `bvartools-dev` skills cover the two layers below.

## 1. Decide where the change belongs

| The change is... | It goes to |
| --- | --- |
| A new sampler, prior inside a sampler, or likelihood | **BayesTS**, then vendored into bvartools (see `bvartools-dev`), then used here. Never a sampler written in R. |
| Heavy numerical work over draws (posterior simulation, a filter, drawing from a closed form) | **BayesTS**. R does the interface and the light tasks; BayesTS works on the model files themselves rather than R shuttling draws model by model. |
| A bug found through bgvars in bvartools behaviour | **Fixed in bvartools**, which is the user's package too, not worked around here. |
| Machinery for working on models stored in HDF5 that is not specific to a *global* model | Written so it can move to bvartools, which is meant to get the same file-first workflow. |
| Anything specific to a global model | Here. |

When a task seems to need a workaround for something that belongs in BayesTS or
bvartools, **propose the upstream change to the user first**. `DESCRIPTION`
pins `bvartools (>= ...)` and declares `Remotes: franzmohr/bvartools`; bump the
minimum when this package starts relying on something new there.

## 2. Anatomy

Classes:

```text
gvarmodel / gvecmodel      the global model in memory: $global, $weights, $submodels
  $submodels[[country]]    a 'modellist' of candidates
    varxsubmodel           inherits 'bvarmodel'
    vecxsubmodel           inherits 'bvecmodel'
gvarexport                 handle of an export folder: $folder, $manifest, $rclass
submodeldata               the input data, one element per country
```

A verb usually has up to four methods:

```text
R/add_foo.gvarmodel.R      maps over every candidate of every sub-model, then restores the class
R/add_foo.gvecmodel.R      the same for a GVEC (often delegating to shared helpers)
R/add_foo.varxsubmodel.R   only when a sub-model differs from a bvarmodel (two blocks of
                           unmodelled variables, foreign-variable naming); otherwise
                           NextMethod() / the bvartools method does the work
R/export_methods.R         the 'gvarexport' method, a one-liner over map_submodels()
```

Mapping over sub-models goes through the helpers in `R/submodel_models.R`
(`.submodel_models()`, `.apply_to_submodel_models()`), which also carry `cores`
and number models for seeds. A method for the export must go through
`map_submodels()`, which reads, applies, writes through a temporary file and
rebuilds the manifest; do not read models into a list.

Global-model invariants that code must keep:

- **Deterministic terms are built once** in `create_gvarmodel()` /
  `create_gvecmodel()` on the global time axis; sub-models take theirs from it,
  so the trend regressor is common to all.
- **Lag orders of unmodelled blocks count from lag zero**, and VECX lag orders
  are those of the model in levels. `vec_to_var.vecxsubmodel()` and
  `submodels_to_gvar()` read coefficient positions in the order *endogenous
  lags, weakly exogenous from lag 0, global from lag 0, deterministic*.
- **Seeds** are numbered by position in the whole manifest, so a run over a
  subset draws what a full run draws. Anything that changes model order breaks
  reproducibility.
- **The export is found through the manifest**, never by walking the directory
  tree.
- **Memory**: a global model is large in its draws. Anything new that touches
  draws should work one model, or one chunk of draws, at a time. Estimating in
  batches with `write_to_hdf5(append = TRUE)` is the standard workflow, not a
  fallback.

Refusals are statistical and explain themselves (SSVS on a VECX sub-model,
structural with an unrestricted covariance, solving from discounted
closed-form sub-models). Error messages name the argument and say what to do.

## 3. Documentation

Roxygen only; `man/*.Rd`, `NAMESPACE` and `Collate` are generated. Check the
installed roxygen2 against `Config/roxygen2/version` in `DESCRIPTION`, then

```text
& "D:/R/R-4.6.1/bin/Rscript.exe" -e "devtools::document()"
```

House style for `@details` is explanatory prose -- what the method does, why,
and what can go wrong -- including a `@section` when a function has a mode worth
its own heading (see `write_to_hdf5.gvarmodel.R`, `submodels_to_gvar.R`).
Several older help pages (`add_priors.varxsubmodel`, `add_priors.vecxsubmodel`)
still describe default values that the code no longer has; fix such text when
touching those files.

## 4. Tests

- testthat 3e, one `tests/testthat/test-<topic>.R` per topic.
- Build objects through the fixtures in `helper-fixtures.R` (`gvar_object()`,
  `gvec_object()`, the estimated fixtures): three countries, two variables,
  tens of draws. Keep new fixtures that small.
- Test invariants and estimands, not only shapes: a stacked global model that is
  wrong still has the right dimensions. Compare against a hand-built solution or
  a known published number where there is one (the GVEC vignette reproduces
  Dees et al. 2007).
- `test-agent-docs.R` runs every ```` ```r ```` block of the user-facing skill
  against the installed package. Use ```` ```text ```` for code that must not run
  (shell commands, paths on this machine).

```text
& "D:/R/R-4.6.1/bin/Rscript.exe" -e "devtools::test(filter = 'export-pipeline')"
& "D:/R/R-4.6.1/bin/Rscript.exe" -e "devtools::test()"
```

`test-agent-docs.R` reads `system.file()`, so it tests the **installed** copy:
run `devtools::install(upgrade = FALSE)` first after editing `inst/agents/`.
`R` is not on `PATH`; call it by full path. Rtools is at `D:/rtools45`, pandoc
inside RStudio at `D:/rstudio`, the BayesTS executable at
`D:/BayesTS/bin/bayests.exe`.

## 5. Vignettes are precompiled

Sources are `vignettes/<name>.Rmd.orig`; the `.Rmd` beside each and the figures
under `vignettes/figures/` are generated and committed. Edit the `.orig`, then

```text
& "D:/R/R-4.6.1/bin/Rscript.exe" vignettes/precompile.R            # all
& "D:/R/R-4.6.1/bin/Rscript.exe" vignettes/precompile.R GVEC       # one
```

The sources deliberately use `devtools::load_all(".")` rather than
`library(bgvars)`; do not propose changing that. Each keeps
`fig.path = "figures/<name>-"`. The GVEC vignette calls BayesTS and writes to
`D:/models/dees2007/`, so it only rebuilds on this machine.

## 6. Finishing a change

A change is done when all of these agree with it:

1. **Roxygen** updated and `devtools::document()` run.
2. **Tests** added or changed, and the relevant files pass.
3. **The user-facing skill** in `inst/agents/` (`AGENTS.md`, `SKILL.md`,
   `references/*.md`) if what a user sees changed: a new function or argument,
   a renamed one, a new refusal, a changed default.
4. **Vignettes** regenerated if their output changed.
5. **The JSS manuscript** at `D:/syncthing/projects/bgvars_jss`
   (`2026-bgvars.tex`, code in `article.R`) when an algorithm or the documented
   workflow changed. A change to an algorithm is not finished until the
   documentation and the paper describe it.
6. **The containerised check** before a push: `D:\workspace-r\r-ci-docker\rcheck.ps1 bgvars`
   (or `./rcheck.sh bgvars` from Git Bash), which mounts the local bvartools and
   should end in `OK`. A host `devtools::check()` does not substitute for it.
   Report its result; do not wire it into a git hook.

There is no `NEWS.md` yet; if one is added, follow the bvartools convention (a
bold sentence saying what the user can now do, then prose).

## 7. Git workflow

Trunk-based: `main` is the development version and may be ahead of the latest
release; releases are tags plus GitHub releases; there is no `dev` branch.

1. Cut a short branch from `main` for one piece of work.
2. Commit with an imperative, plain-English title saying what changed for the
   user (`Solve a global model from its export, a piece of the chain at a
   time`), and a body of prose paragraphs on why and on what was decided
   against. End with the `Co-Authored-By` line the session provides.
3. Merge back with `git merge --no-ff` and a title `Merge: <lower-case summary>`,
   then delete the branch.

Commit only when asked. Pushing is the user's; do not remind them about unpushed
commits or how far `main` is ahead of `origin`.
