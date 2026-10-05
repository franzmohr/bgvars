# Specifying a global model

Everything before the posterior simulation: the data, the weights, the
sub-models and their priors. Every `r` block here runs in the package's test
suite.

## The data: `submodeldata`

A global model starts from a named list of class `submodeldata`, one element per
country, each holding

- `endogen` — a `ts` of the country's own series, variables in columns;
- `weights` — a `ts` (usually annual) of the country's trade with every other
  country, one column per country.

Global series (the oil price, commodity prices) are a separate multivariate `ts`.
The shipped data sets come in that shape:

| Data set | Content |
| --- | --- |
| `gvar2023` | Mohaddes and Raissi (2024): 33 countries, 1979Q2–2023Q3, `y`, `Dp`, `r`, `lr`, `ep`, `eq`; global `poil`, `pmat`, `pmetal`; `region_weights` (PPP-GDP) |
| `gvar2019` | The 2020 vintage of the same database, to 2019Q4 |
| `dees2007` | Dees, di Mauro, Pesaran and Smith (2007): 26 economies with the euro area as one, 1979Q2–2003Q4, `y`, `Dp`, `rs`, `rl`, `q`, `ep`; global `poil` |

Not every country has every variable. A variable listed in `endogen` or `exogen`
that a country lacks is simply left out of that country's sub-model.

```r
library(bgvars)
data("dees2007")
submodel_data <- dees2007[["submodel_data"]]

stopifnot(inherits(submodel_data, "submodeldata"),
          all(c("endogen", "weights") %in% names(submodel_data[["US"]])))

# Keep a few countries; the class survives
small <- select_list_elements(submodel_data, c("EA", "JP", "US"))
stopifnot(inherits(small, "submodeldata"), length(small) == 3)

# Differences of chosen variables, times 100
growth <- diff(small, variables = "y", multi = 100)
```

`create_regions(submodel_data, region_weights, regions, period)` aggregates
countries into regions (e.g. a euro area of eight members, weighted by PPP-GDP)
and returns a new `submodeldata`. Whether to form regions is a modelling choice;
every country can equally be its own sub-model.

## Weights

`add_weight_matrices(object, submodel_data, period)` builds, per sub-model, the
matrix `W_i` that maps the vector of all endogenous variables of the global model
to the sub-model's own and foreign variables.

| `period` | Weights |
| --- | --- |
| a vector of periods, e.g. `1999:2001` | constant, from the trade summed over those periods |
| one integer that is not a period of the weight series, e.g. `3` | time varying, rolling sums over the last 3 periods |
| one period of the series, e.g. `2001` | constant, from that period alone |

Periods the weight series does not contain are refused.
`get_weight_matrix(object, "US")` returns a sub-model's matrix.

## Global model object

`create_gvarmodel(submodel_data, global_data)` and
`create_gvecmodel(submodel_data, global_data)` differ only in the class they
return, which decides what `add_submodels()` builds (VARX or VECX). The object
has three elements: `global` (the data of every country, the global series, the
deterministic terms and an index), `weights` and `submodels`. The deterministic
terms are built once on the global time axis, so every sub-model shares the same
trend regressor whatever its lag order.

## Sub-models

`add_submodels(object, ...)` builds a sub-model for **every** country with the
same arguments. For country-specific specifications, build them one at a time
with `create_varxsubmodel()` or `create_vecxsubmodel()` and assign the result:

```r
object <- create_gvecmodel(submodel_data = small,
                           global_data = dees2007[["global_data"]])
object <- add_weight_matrices(object, submodel_data = small, period = 1999:2001)

for (country in names(small)) {
  object[["submodels"]][[country]] <- create_vecxsubmodel(object,
    submodel = country,
    endogen = c("y", "Dp", "rs"), p_endogen = 2,
    exogen = c("y", "Dp", "rs"), p_exogen = 2,
    r = if (country == "US") 2 else 1,
    const = "unrestricted", trend = "restricted",
    iterations = 100, burnin = 50)
}

stopifnot(inherits(object[["submodels"]][["US"]], "modellist"),
          inherits(object[["submodels"]][["US"]][[1]], "vecxsubmodel"))
```

What a sub-model is:

- `object$submodels[[country]]` is always a `modellist`, even with one
  candidate. `object$submodels[[country]][[1]]` is the model.
- Foreign (weakly exogenous) variables carry the suffix `.s`: `y.s`, `Dp.s`.
  Regressor names add the lag, `y.s.l00` for the contemporaneous term.
- `get_submodel_specifications(object$submodels[[country]])` tabulates the
  candidates: type, variables and lag orders of each block, rank, variable
  selection.

```r
spec <- get_submodel_specifications(object[["submodels"]][["US"]])
stopifnot(spec$type == "VECX", object[["submodels"]][["US"]][[1]][["model"]][["rank"]] == 2)
```

### Lag orders

| Argument | VARX sub-model | VECX sub-model |
| --- | --- | --- |
| `p_endogen` | lags of the own variables | lag order of the model in levels; `p_endogen - 1` lagged differences |
| `p_exogen` | blocks of foreign variables from lag 0: `1` = contemporaneous only, `2` adds lag 1, `0` = none | the same count in levels, at least 1; the lagged level always enters the error correction term, and `p_exogen - 1` differences from `Δx_t` on |
| `s` | the same count for the global variables; `NULL` = none | the same, at least 1 |

A VARX\*(p, q) of the GVAR literature is `p_endogen = p, p_exogen = q + 1` as a
VECX sub-model. The literature searches p and q over 1–2.

### Deterministic terms

- VARX: `deterministic = "none"`, `"const"` (default), `"trend"`, `"both"`, and
  `seasonal = TRUE`.
- VECX: `const`, `trend` and `seasonal`, each `"restricted"` (in the
  cointegration space), `"unrestricted"` or `NULL`. Case IV of Pesaran, Shin and
  Smith (2000), used by Dees et al., is `const = "unrestricted",
  trend = "restricted"`.

### Global variables and the dominant unit

`global = "poil", s = 1` puts the oil price into every sub-model as an
exogenous regressor. In the Dees et al. (2007) setup the oil price is instead
**endogenous in the US model**: add it to the US `endogen` series before
creating the model, list it among the US model's `endogen`, and list `poil`
among the `exogen` of every other country, which then takes it as a foreign
variable. `references/recipes.md` shows the code.

### Model types

The arguments are bvartools' and mean the same:

| Argument | Values |
| --- | --- |
| `error` | `"wishart"`, `"gamma"`, `"gamma+covar"`, `"sv"`, `"sv+covar"` |
| `tvp` | `TRUE` for random-walk coefficients |
| `varsel` | `"none"`, `"bvs"`, `"ssvs"` (VARX only) |
| `structural` | `TRUE` for contemporaneous coefficients among the own variables |
| `algorithm` | `NULL`; `"KLGS2010"` (VECX); `"discount"` for the closed-form discounted model, with `delta_beta`, `delta_sigma` in `(0, 1]` |

`align_model_obs(object)` puts every candidate of every sub-model on the sample
they share, which comparisons of lag orders need.

## Priors

`add_priors(object, coef = ..., sigma = ..., coint = ..., varsel = ...)` passes
the lists to every sub-model. Nothing has a default.

| Sub-model | `coef` | `sigma` | `coint` |
| --- | --- | --- | --- |
| constant, Wishart | `v_i` (or `minnesota`), optional `v_i_det`, `const` | `df`, `scale` | VECX: `v_i`, `p_tau_i` |
| `error = "gamma"` | as above | `shape`, `rate` | as above |
| `error = "sv"`, `"sv+covar"` | as above | `shape`, `rate`, `mu`, `v_i`, `state_variance`, `offset` | as above, optionally `g_i` for the loadings' prior |
| `tvp = TRUE` | add `shape`, `rate`, optional `rate_det`; VECX also `rate_alpha` for the loadings | as above | VECX: `rho` for time varying cointegration |
| `varsel = "bvs"` | as above | as above | pass `varsel = list(inprior = 0.5, exclude_det = TRUE)` |

`varsel$covar` says whether selection reaches the error covariance block, and
every sub-model with such a block answers it, whether or not its coefficients
are time varying. Left out it is `FALSE`, so the coefficients are selected over
and the covariances are not, and a run that wants selection there asks for it.
Either way it is worth a thought: a generalised impulse response is built from
the off-diagonals of the error covariance, so selecting them to zero collapses
it towards an orthogonalised response, and the block carries only
`k_i(k_i-1)/2` elements per sub-model against the coefficients' many more.

```r
object <- add_priors(object,
                     coef = list(v_i = 0),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 3, scale = 1e-4))

bad <- tryCatch(add_priors(object, coef = list(v_i = 0), sigma = list(df = 3, scale = 1e-4)),
                error = function(e) "refused")
stopifnot(identical(bad, "refused"))   # a GVEC needs 'coint'
```

A flat coefficient prior (`v_i = 0`), a flat cointegration space prior
(`coint = list(v_i = 0, p_tau_i = 1)`) and `scale_error_correction()` are what
constant-coefficient GVEC models are estimated with in the vignettes. For time
varying sub-models the random-walk rates are the decisive choice: a `rate` of
`1e-5` or smaller for the coefficients, and a `rate_alpha` several orders of
magnitude smaller still, keeps the drift from absorbing the residuals (see the
TVP-SV vignettes).

`minnesota_prior()` and `inclusion_prior()` have methods for single sub-models
when the prior should be inspected before it is used.

## What is refused, and why

- **SSVS on a VECX sub-model** — not implemented for error correction
  sub-models; use `varsel = "bvs"`. The loadings are never subject to
  selection, so the error correction term stays in the model.
- **A structural sub-model with `error = "wishart"`, `"gamma+covar"` or
  `"sv+covar"`** — `A_0` is not identified against an unrestricted covariance;
  use `"gamma"` or `"sv"`.
- **`algorithm = "discount"`** needs `error = "wishart"`, `burnin = 0`,
  `thin = 1`, no variable selection and no structural model. It is estimated by
  BayesTS on its file and by nothing in R, and its closed-form posterior has no
  draws, so `submodels_to_gvar()` refuses it: it is for comparing
  specifications by their exact marginal likelihood, not for the model the
  impulse responses come from.
- **`submodels_to_gvar()` with several candidates for a sub-model** — choose one
  per country first.
