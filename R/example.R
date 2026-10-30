spt_example_data <- function(n_source = 600L, n_target = 300L, seed = 1L) {
  if (any(!is.finite(c(n_source, n_target))) || any(c(n_source, n_target) < 20) ||
      any(c(n_source, n_target) != as.integer(c(n_source, n_target))))
    .stop("Example cohort sizes must be integers of at least 20.")
  .with_seed(seed, {
    generate <- function(n, target) {
      age <- stats::runif(n, 40 + 5 * target, 80 + 5 * target)
      sex <- stats::rbinom(n, 1, .5); biomarker <- stats::rnorm(n, .2 * target)
      mu <- 3.4 - .04 * (age - 60) - .25 * sex + .2 * biomarker^2 + .15 * target
      event_time <- exp(mu + stats::rlogis(n, scale = .7))
      censor_time <- stats::rexp(n, rate = .014)
      data.frame(time = pmin(event_time, censor_time), status = as.integer(event_time <= censor_time),
                 age = age, sex = sex, biomarker = biomarker)
    }
    list(source = generate(n_source, 0), target = generate(n_target, 1))
  })
}
