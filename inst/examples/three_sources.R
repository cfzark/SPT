# Run after installing SPT. All cohorts are synthetic.
library(SPT)
d <- spt_example_data(600, 300, seed = 42)
target_formula <- ~ sex + s(age, bs="ts", k=5) + s(biomarker, bs="ts", k=5)
sources <- list(
  GAM = spt_source(d$source, "time", "status",
    ~ sex + s(age, bs="ts", k=10) + s(biomarker, bs="ts", k=10),
    learner="gam", horizon=24),
  Cox = spt_source(d$source, "time", "status", ~ sex+age+biomarker,
    learner="cox"),
  NN = spt_source(d$source, "time", "status", ~ sex+age+biomarker,
    learner="nn", horizon=24)
)
models <- lapply(sources, function(source) {
  spt_fit(d$target, "time", "status", 24, target_formula, source)
})
print(lapply(models, function(m) list(rho=m$rho, tau=m$tau, cv=m$cv$table)))
print(sapply(models, predict, newdata=head(d$target)))
