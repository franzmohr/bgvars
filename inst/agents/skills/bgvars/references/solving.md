# Solving and analysing the global model

Every `r` block here runs in the package's test suite.

## What `submodels_to_gvar()` does

Sub-model `i` is a VARX\* model in its own variables `x_it` and its foreign
variables `x*_it = W_i`-weighted averages of the other countries' variables.
Written in the vector `y_t` of every country's endogenous variables, its
contemporaneous terms are `(A_0i, -A*_0i) W_i y_t`; stacking them over all
sub-models gives the matrix `G`, and the lagged terms stack the same way:

    G y_t = A_1 y_{t-1} + ... + A_p y_{t-p} + C d_t + D g_t + u_t

The function solves this **draw by draw** for the reduced form
`y_t = G^{-1} A_1 y_{t-1} + ... + G^{-1} u_t` and returns it as a `bvarmodel`:

| Element | Holds |
| --- | --- |
| `model$k`, `model$endogen` | the number and names of all endogenous variables, `<country>_<variable>` |
| `model$p` | the largest lag order of any sub-model |
| `model$algorithm` | `"VarNormalWishart"` — whatever the sub-models were, the solved model has constant coefficients and a full error precision, which is what bvartools dispatches on |
| `posterior$a$coeffs`, `posterior$u_sigma_inv$coeffs` | draws of the reduced form, draws in rows, as for any bvartools model |
| `g` | the draws of `G`, one row per draw, column-major `k x k` |
| `data$original$exogen` | the global series, to be continued for a forecast |

Properties that follow:

- **Exactly one model per sub-model**; with candidates left the function stops.
- **Sub-models may differ** in lag orders, deterministic terms and even type: a
  VECX sub-model is first rewritten in levels by `vec_to_var()`, so a global
  model mixing VARX and VECX sub-models solves in one step. The rank plays no
  further role.
- **The structural error covariance is block diagonal**, one block per
  sub-model, so the reduced-form errors `G^{-1} u_t` are correlated across
  countries through `G` only.
- **`period`** picks the date at which time varying weights, TVP coefficients
  and SV covariances are taken; default the last period of the sample. It is a
  **row number of the global data**, not a time. For a path over time, call the
  function once per period.
- A **scaled** GVEC (`scale_error_correction()`) must be rescaled first.
- The draws are not checked for stability. A drawn global model can have roots
  outside the unit circle, especially in GVECs with diffuse cointegration
  priors; check the companion matrix before summarising impulse responses.

```r
library(bgvars)
set.seed(7)
data("gvar2023")
submodel_data <- select_list_elements(gvar2023[["submodel_data"]], c("DE", "JP", "US"))
object <- create_gvarmodel(submodel_data, global_data = gvar2023[["global_data"]])
object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
object <- add_submodels(object, endogen = c("y", "Dp", "r"), p_endogen = 2,
                        exogen = c("y", "Dp", "r"), p_exogen = 2,
                        global = "poil", s = 1,
                        iterations = 200, burnin = 100)
object <- add_priors(object, coef = list(v_i = 0), sigma = list(df = 3, scale = 1e-4))
object <- add_initial_values(object)
object <- add_posterior_coefficients(object)
gvar <- submodels_to_gvar(object)

k <- gvar[["model"]][["k"]]
p <- gvar[["model"]][["p"]]
draws <- nrow(gvar[["posterior"]][["a"]][["coeffs"]])
stopifnot(k == 9, p == 2, draws == 200,
          dim(gvar[["g"]]) == c(draws, k * k),
          "US_r" %in% gvar[["model"]][["endogen"]])

# Largest modulus of the companion matrix, draw by draw
max_modulus <- apply(gvar[["posterior"]][["a"]][["coeffs"]][, 1:(k * k * p)], 1, function(a) {
  a <- matrix(a, k)
  companion <- rbind(a, cbind(diag(k * (p - 1)), matrix(0, k * (p - 1), k)))
  max(Mod(eigen(companion, only.values = TRUE)$values))
})
stable <- max_modulus < 1
```

## Impulse responses and variance decompositions

The solved model is a `bvarmodel`, so bvartools' `irf()` and `fevd()` apply. In
a GVAR there is no natural recursive ordering across countries, so the usual
choice is the **generalised** impulse response of Pesaran and Shin (1998):

```r
gir <- irf(gvar, impulse = "US_r", response = "DE_y", n_ahead = 12,
           ci = 0.68, type = "gir", shock = "sd")
gfevd <- fevd(gvar, response = "DE_y", n_ahead = 12, type = "gir",
              normalise_gir = TRUE)

stopifnot(inherits(gir, "bvarirf"), nrow(gir) == 13)
```

`shock = "sd"` scales the shock to one standard deviation of the impulse
variable's reduced-form residual. **What is shocked is the reduced-form error
`G^{-1} u_t`**, with the covariance the solution assembles from the block
diagonal country covariances. The generalised responses of Dees et al. (2007,
eq. 21) shock the stacked country residuals `u_t` with their full covariance,
cross-country blocks included; that is a different object, and the GVEC
vignette shows how to compute it from `gvar$g` and the residuals.

## Forecasts

A forecast of the global model needs the future path of the global variables:

```r
poil <- gvar[["data"]][["original"]][["exogen"]]
tsp_poil <- stats::tsp(poil)
future <- stats::ts(matrix(rep(poil[nrow(poil), ], each = 5), 5, ncol(poil)),
                    start = tsp_poil[2], frequency = tsp_poil[3])
dimnames(future)[[2]] <- dimnames(poil)[[2]]

gvar <- add_forecast_input(gvar, n_ahead = 4, exogen = future)
gvar <- add_posterior_forecasts(gvar)
fcst <- predict(gvar, n_ahead = 4)
```

Holding the oil price at its last value is only a placeholder; supply the path
the analysis assumes.

## Working on one period of a time varying model

```text
global_time <- stats::time(object[["global"]][["endogen"]])
gvar_2008q4 <- submodels_to_gvar(object, period = which(abs(global_time - 2008.75) < 1e-6))
```

`period` is the position of a period in the global data (`1` is its first
quarter) and must lie inside the estimation sample of the time varying
sub-models. Constant sub-models contribute the same draws in
every period.
