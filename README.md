# SPT: Source Projection Transfer

SPT transfers an external survival-risk score into a user-specified target GAM.
This package implements the calibrated **NoResidual/SPT** estimator with a
profiled shape-transfer weight `rho` in `[0,1]` and fixed or cross-validated
soft-label weight `tau`. It does not implement DST, raw-teacher variants,
residual-channel models, or the paper's comparison methods.

## Installation

```r
install.packages(c("mgcv", "survival", "nnet"))
install.packages("SPT_0.1.0.tar.gz", repos = NULL, type = "source")
library(SPT)
```

Alternatively, `R CMD INSTALL packages/SPT` from the development project.
No compiler, Slurm installation, private data, or project-specific paths are needed.
This is a local development release, not a CRAN release.

## GAM source and user-specified target

```r
d <- spt_example_data()
source <- spt_source(
  d$source, time = "time", status = "status", horizon = 24,
  learner = "gam",
  formula = ~ sex + s(age, bs = "ts", k = 10) +
                   s(biomarker, bs = "ts", k = 10)
)
target_formula <- ~ sex + s(age, bs = "ts", k = 5) +
                         s(biomarker, bs = "ts", k = 5)
fit <- spt_fit(d$target, "time", "status", horizon = 24,
               target_formula = target_formula, source = source)
print(fit)
predict(fit, head(d$target))
fit$cv$table
fit$rho
```

Default `tau=NULL` uses `{0, .1, .5, 1}` and three-fold **paired one-SE-half**.
The reference `.5` is retained unless the empirical winner improves target
validation log loss by more than one paired SE. Set `tau=.5` for a fixed fit,
or `selection="min_cv"` for minimum validation loss. User grids may exceed one.
`tau=0` still borrows through the shape penalty; it is not target-only GAM.

## Cox and NN sources

```r
cox_source <- spt_source(d$source, "time", "status",
                         ~ age + sex + biomarker, learner = "cox")
cox_spt <- spt_fit(d$target, "time", "status", 24,
                   target_formula, cox_source)
# Reuse this Cox source at 30/36: absolute risk is transformed to a logit.

nn_source <- spt_source(d$source, "time", "status",
                        ~ age + sex + biomarker, learner = "nn", horizon = 24)
nn_spt <- spt_fit(d$target, "time", "status", 24,
                  target_formula, nn_source)
```

NN defaults: one hidden layer, size 4/8, decay .1/1/10, five-fold source-only
CV, and three deterministic starts averaged in probability space. GAM and NN
sources need a separate fit for each horizon. Cox needs one fit. Source GAM
smooth complexity is entirely controlled by its formula; `k=10` is an example,
not a hidden override. NN has no spline `k`.

## Target structure belongs to the user

```r
interaction_formula <- ~ sex + s(age, bs="ts", k=5) +
  s(biomarker, bs="ts", k=5) +
  ti(age, biomarker, bs=c("ts", "ts"), k=c(3,3))
```

The package neither chooses interactions nor removes variables. Choose the
formula on training data; if structure selection uses target outcomes, perform
it within an appropriate nested validation scheme outside this API. Formula
inputs are one-sided, with an intercept and explicitly named variables.
Ordinary smooths, tensor interactions, and factor-by smooths are supported;
linked `id` smooths and offsets are not supported in v0.1.0. All variables must
be complete, factors must use compatible levels, and the target basis must
have full column rank. Source predictors must be available wherever the source
is evaluated on target training rows. Deployment needs only target predictors.

Continuous raw predictors are standardized within each fitting sample before
formula evaluation. Thus formula transformations refer to standardized inputs;
precompute a transformed variable if raw-unit transformation semantics matter.

## What is fitted

1. A training-only IPCW quasibinomial target GAM, fitted by REML, supplies `B`
   and the fixed smoothness penalty `S`.
2. A Gaussian REML GAM projects the source logit into the target function space.
   Its fitted values are mapped to coefficients `c` in the target design. This
   is a **smooth projection**, not an exact orthogonal least-squares projection.
3. An IPCW intercept calibration of the original source logit supplies `q`.
4. Optimize target Bernoulli loss + `tau` times soft-label loss +
   `(theta-rho*c)' S (theta-rho*c)/2`, profiling `rho` analytically in `[0,1]`.
5. For automatic tau, repeat all target nuisance estimation inside each fold,
   select by paired validation loss, then refit on the full target training set.

`predict(fit, ..., type="components")` returns logit-scale target-only,
projection, target-correction, borrowed-projection and SPT predictions. These
include intercepts; they are not centered partial effects. The last two SPT
channels sum to the SPT logit. No source object is required at prediction time.

## Numerical behavior and reproducibility

- Historical defaults: REML/gamma 1.4, training IPCW cap 20, censor-survival
  floor .001 and probability clipping .001. Censor KM is right-continuous,
  including ties. Validation uses training-only KM without capping weights.
- Independent censoring is assumed; no conditional censoring model is provided.
- Failed GAM fits or tau candidates stop with an error. There is no automatic
  formula reduction, scale change, or silent deletion of a failed candidate.
- The SPT optimizer can retry the same objective with `nlminb`, reporting a
  warning and retaining diagnostics. KKT and objective checks remain enforced.
- `spt_control(gam_scale=1)` explicitly changes the quasibinomial nuisance fit.
  It is an opt-in alternative, not an automatic reproduction of default REML.
- Set `fold_id` for exact user partitions. Seeds do not alter the caller RNG.
- Separate fixed-horizon fits do not guarantee monotone risks over time.
- CV uncertainty is a selection heuristic, not a patient-level confidence bound.

The package's generic fold generator is not the historical clinical partition
generator. Reproducing a paper experiment also requires its exact folds,
formulas, observations, and preprocessing; those belong in the separate paper
replication archive. All bundled examples/tests use synthetic data.

Fitted `mgcv`, Cox and NN objects, model frames and CV losses may contain
training information. Do not publish real-data model RDS files without review.

## Development checks

```sh
R CMD build packages/SPT
R CMD check --no-manual SPT_0.1.0.tar.gz
```

Author/maintainer: Fangzhou Cao <caofzh@mail3.sysu.edu.cn>.
License: GPL-3. Publication name: Source Projection Transfer (SPT).
