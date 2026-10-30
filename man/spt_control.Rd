\name{spt_control}
\alias{spt_control}
\title{Numerical and IPCW Controls}
\usage{spt_control(gamma = 1.4, eps_g = .001, weight_cap = 20,
  prob_clip = .001, kkt_tol = 2e-6,
  gam_scale = 0, nn_maxit = 2000L)}
\arguments{
\item{gamma}{GAM smoothing inflation parameter.}
\item{eps_g}{Lower bound for censoring survival probabilities.}
\item{weight_cap}{Cap on training IPCW weights; \code{Inf} disables capping.}
\item{prob_clip}{Prediction and soft-label probability clipping threshold.}
\item{kkt_tol}{Maximum accepted optimizer stationarity residual in whitened coordinates.}
\item{gam_scale}{Quasibinomial GAM scale: zero estimates scale; a positive value
fixes it. This changes the nuisance estimator and is never activated silently.}
\item{nn_maxit}{Maximum NN optimizer iterations.}}
\description{Construct controls shared by source and target fits.}
\details{The historical IPCW convention uses the right-continuous censoring
Kaplan--Meier estimate at min(observed time, horizon), including ties. Validation
uses a training-only censoring KM and the probability floor, without weight
capping. Independent censoring is assumed; conditional censoring models are
not implemented. GAM fits use REML. Failed fits stop rather than silently
changing scale, removing terms, or dropping candidates.}
\value{A named list.}
