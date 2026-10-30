\name{spt_source}
\alias{spt_source}
\alias{predict.spt_source}
\alias{print.spt_source}
\title{Fit a GAM, Cox, or IPCW Neural-Network Source Learner}
\usage{
spt_source(data, time, status, formula,
  learner = c("gam", "cox", "nn"),
  horizon = NULL, folds = 5L, fold_id = NULL, seed = 20261029L,
  nn_grid = expand.grid(size = c(4L, 8L), decay = c(.1, 1, 10)),
  nn_starts = 3L, control = spt_control())
\method{predict}{spt_source}(object, newdata, type = c("response", "link"),
  horizon = object$horizon, ...)
\method{print}{spt_source}(x, ...)
}
\arguments{
\item{data}{Source data frame with complete predictors and right-censored outcomes.}
\item{time, status}{Column names for observed follow-up and event indicator (1 = event).}
\item{formula}{One-sided predictor formula. GAM accepts mgcv smooth terms;
Cox and NN accept ordinary model-matrix terms. No outcome columns or dot expansion.}
\item{learner}{Source learner type.}
\item{horizon}{Positive prediction time, in the units of \code{time}. Required for
GAM and NN. A Cox fit can be reused at different horizons within observed follow-up.}
\item{folds, fold_id}{Number of source-only NN tuning folds, or explicit fold labels.}
\item{seed}{Local RNG seed; the calling RNG state is restored.}
\item{nn_grid}{Data frame with positive integer \code{size} and nonnegative
\code{decay}. Selected by source-only IPCW validation log loss.}
\item{nn_starts}{Number of fitted networks whose probabilities are averaged.}
\item{control}{Controls from \code{\link{spt_control}}.}
\item{object, x}{A fitted source object.}
\item{newdata}{Prediction covariates with matching names and factor levels.}
\item{type}{Raw risk probability or the logit score used by SPT.}
\item{...}{Reserved for S3 compatibility.}}
\description{Train an external source model independently of target outcomes.}
\details{GAM fits IPCW quasibinomial smooths using REML. Numeric predictors
with more than two distinct values are standardized on the fitting sample
before formula evaluation (including transformations). For GAM the link is
the raw fitted logit. Cox uses Efron partial likelihood and the estimated
baseline cumulative hazard to form absolute risks: its link is the clipped
absolute-risk logit, not the Cox linear predictor. NN uses one sigmoid hidden
layer, weighted Bernoulli loss, weight decay, source-only CV, and a deterministic
multi-start probability ensemble. Its design columns are standardized within
each training fold; its link is the clipped probability logit. Source and target
samples must be disjoint. Complete-case selection or imputation is the user's
responsibility. New factor levels are not supported. Objects contain fitted
training information and are not automatically safe to publish.}
\value{An \code{spt_source} object containing the fitted learner, scaling,
horizon, and (for NN) tuning results. Prediction returns a numeric vector.}
\examples{
d <- spt_example_data(300, 150)
src <- spt_source(d$source, "time", "status",
  ~ sex + s(age, bs = "ts", k = 5) + s(biomarker, bs = "ts", k = 5),
  learner = "gam", horizon = 24)
predict(src, head(d$target))
cx <- spt_source(d$source, "time", "status", ~ age + sex + biomarker,
                 learner = "cox")
predict(cx, head(d$target), horizon = 30)
}
