\name{spt_fit}
\alias{spt_fit}
\alias{predict.spt}
\alias{coef.spt}
\alias{print.spt}
\title{Fit Calibrated Source Projection Transfer}
\usage{
spt_fit(data, time, status, horizon, target_formula, source,
  tau = NULL, tau_grid = c(0, .1, .5, 1),
  selection = c("one_se_half", "min_cv"), reference_tau = .5,
  folds = 3L, fold_id = NULL, seed = 20261029L,
  control = spt_control())
\method{predict}{spt}(object, newdata,
  type = c("response", "link", "components"), ...)
\method{coef}{spt}(object, ...)
\method{print}{spt}(x, ...)
}
\arguments{
\item{data}{Target training data frame with complete predictors.}
\item{time, status}{Names of observed-time and event (0/1) columns.}
\item{horizon}{One positive fixed prediction horizon.}
\item{target_formula}{User-chosen one-sided mgcv formula, including an intercept.
For example a linear sex term and \code{s(age, bs="ts", k=5)}. Interactions such as
\code{ti(x1,x2)} or factor-by smooths can be specified explicitly. No automatic
structure selection, offsets, dot expansion, or linked smooth IDs.}
\item{source}{Independent fitted \code{\link{spt_source}}.}
\item{tau}{Fixed nonnegative soft-label strength. \code{NULL} selects by target CV.}
\item{tau_grid}{Unique nonnegative candidate strengths.}
\item{selection}{Paired one-SE reference rule or empirical minimum CV loss.}
\item{reference_tau}{Preferred value for the paired one-SE rule; must occur in
the grid. Default 0.5.}
\item{folds, fold_id}{Number of target CV folds or explicit fold labels.}
\item{seed}{Local fold-generation seed.}
\item{control}{Controls from \code{\link{spt_control}}.}
\item{object, x}{Fitted SPT object.}
\item{newdata}{Target covariates for prediction.}
\item{type}{Risk probability, logit score, or logit-scale component predictions.}
\item{...}{Reserved for S3 compatibility.}}
\description{Estimate a source-centered smoothness penalty and calibrated
soft-label transfer model in a user-specified target GAM basis.}
\details{
The target IPCW quasibinomial GAM supplies a design B and a smoothing penalty
S (REML smoothing parameters fixed thereafter). An unweighted Gaussian REML
GAM maps the source logit onto the same target function space; its fitted
values are represented as B c. This is a smooth projection, not an exact
unpenalized orthogonal projection. An IPCW intercept calibration of the source
logit supplies probabilities q.

For fixed tau, SPT minimizes
\deqn{\sum_i w_i\{b(\eta_i)-y_i\eta_i\}+
\tau\sum_i\{b(\eta_i)-q_i\eta_i\}+
\tfrac12(\theta-\rho c)^T S(\theta-\rho c),}
where eta = B theta and b is the Bernoulli log-partition function.
The bounded coefficient rho in [0,1] is profiled analytically as the clipped
ratio \eqn{c^T S\theta/(c^T Sc)}. If the source direction has zero penalty
curvature, rho is defined as zero. This avoids adding a collinear source
projection as an extra target design column. No residual score enters eta.
Tau=0 still permits source-centered shape transfer and is not target-only GAM.

Automatic selection refits the target GAM, projection, calibration, and tau
candidates independently inside each target training fold. Losses use the
training-fold censoring KM. The one-SE-half rule retains the reference unless
its mean paired loss disadvantage relative to the empirical winner exceeds
one paired standard error. Otherwise it selects the empirical winner; ties
in minimum loss favor smaller tau. This is a selection heuristic, not a
confidence guarantee. The chosen tau is refitted on all target training rows.
Standardization is training-only and performed before formula evaluation.
Structure must be fixed by the caller; choosing a formula with held-out
outcomes would leak information. A source fit is held fixed through target CV.
}
\value{An \code{spt} object with coefficients \code{theta}, correction
\code{beta}, \code{rho}, \code{tau}, calibration, target GAM, projection,
penalty matrix, CV losses and fold IDs, and convergence diagnostics.
Prediction with \code{type="components"} returns target-only logits, raw
projected logits, target correction, rho-weighted projection, and SPT logits;
the last two channels sum to the SPT logit. These include intercepts and are
not centered single-variable effects. Prediction needs no source model.
Response predictions are clipped. Fitted objects retain training information.}
\examples{
d <- spt_example_data(400, 220, seed = 3)
src <- spt_source(d$source, "time", "status",
  ~ sex + s(age, bs="ts", k=10) + s(biomarker, bs="ts", k=10),
  horizon = 24)
fit <- spt_fit(d$target, "time", "status", 24,
  ~ sex + s(age, bs="ts", k=5) + s(biomarker, bs="ts", k=5),
  source = src, tau = .5)
print(fit)
predict(fit, head(d$target))
}
