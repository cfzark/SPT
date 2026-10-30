\name{SPT-package}
\alias{SPT-package}
\alias{SPT}
\docType{package}
\title{Source Projection Transfer for Survival Prediction}
\description{Calibrated source projection transfer (SPT) for a single fixed
survival horizon. Fit an external learner with \code{\link{spt_source}}, then
fit a user-specified target GAM using \code{\link{spt_fit}}.}
\details{The package implements SPT (the calibrated NoResidual estimator),
not the historical residual-channel or direct-score transfer ablations.
The user determines target structure. Separate horizon-specific models need
not yield monotone risks across horizons.}
\author{Fangzhou Cao \email{caofzh@mail3.sysu.edu.cn}}
