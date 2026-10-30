\name{spt_example_data}
\alias{spt_example_data}
\title{Generate Small Synthetic Survival Cohorts}
\usage{spt_example_data(n_source = 600L, n_target = 300L, seed = 1L)}
\arguments{\item{n_source, n_target}{Cohort sizes, integers at least 20.}
\item{seed}{Local seed; caller RNG state is restored.}}
\description{Generate independent source and target cohorts with age, sex,
biomarker, log-logistic survival times, and independent exponential censoring.
These are instructional synthetic data, not the manuscript DGP or patient data.}
\value{A list containing \code{source} and \code{target} data frames.}
