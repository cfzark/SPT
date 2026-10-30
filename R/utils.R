utils::globalVariables(".spt_weight")
.stop <- function(...) stop(..., call. = FALSE)
.clip <- function(x, eps) pmin(pmax(x, eps), 1 - eps)
.softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

spt_control <- function(gamma = 1.4, eps_g = .001, weight_cap = 20,
                        prob_clip = .001, kkt_tol = 2e-6,
                        gam_scale = 0, nn_maxit = 2000L) {
  z <- list(gamma = gamma, eps_g = eps_g, weight_cap = weight_cap,
            prob_clip = prob_clip, kkt_tol = kkt_tol,
            gam_scale = gam_scale, nn_maxit = nn_maxit)
  .control(z)
}
.control <- function(x) {
  defaults <- list(gamma = 1.4, eps_g = .001, weight_cap = 20,
                   prob_clip = .001, kkt_tol = 2e-6,
                   gam_scale = 0, nn_maxit = 2000L)
  if (!is.list(x) || length(setdiff(names(x), names(defaults))))
    .stop("Unknown control entry.")
  z <- utils::modifyList(defaults, x)
  if (any(!vapply(z, function(a) is.numeric(a) && length(a) == 1 && !is.na(a), logical(1))))
    .stop("Control values must be numeric scalars.")
  if (!is.finite(z$gamma) || z$gamma <= 0 || z$eps_g <= 0 || z$eps_g >= 1 ||
      z$prob_clip <= 0 || z$prob_clip >= .5 || z$weight_cap <= 0 ||
      !is.finite(z$kkt_tol) || z$kkt_tol <= 0 || !is.finite(z$gam_scale) ||
      z$gam_scale < 0 || !is.finite(z$nn_maxit) || z$nn_maxit < 1)
    .stop("Invalid control values.")
  z
}
.data_check <- function(data, time, status, horizon = NULL) {
  if (!is.data.frame(data) || nrow(data) < 3 || anyDuplicated(names(data)))
    .stop("data must be a data frame with at least three rows and unique names.")
  if (length(time) != 1 || length(status) != 1 || identical(time, status) ||
      !all(c(time, status) %in% names(data))) .stop("Specify distinct time/status columns.")
  y <- data[[time]]; d <- data[[status]]
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0) ||
      !(is.numeric(d) || is.logical(d)) ||
      anyNA(d) || !all(d %in% c(0, 1)))
    .stop("Times must be finite and nonnegative; status must be 0/1 (1 = event).")
  if (!is.null(horizon) && (length(horizon) != 1 || !is.finite(horizon) || horizon <= 0))
    .stop("horizon must be one finite positive number in the same units as time.")
  invisible(TRUE)
}
.formula <- function(formula, data, response = ".spt_response") {
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .stop("Use a one-sided predictor formula, e.g. ~ sex + s(age, bs='ts', k=5).")
  vars <- all.vars(formula)
  if (!length(vars) || !all(vars %in% names(data)) || "." %in% vars)
    .stop("Name predictor columns explicitly; formulas using '.' are not supported.")
  if (any(startsWith(vars, ".spt_"))) .stop("Predictor names beginning with '.spt_' are reserved.")
  if (anyNA(data[vars]) || any(vapply(data[vars], is.character, logical(1))))
    .stop("Predictors must be complete; convert categorical strings to factors first.")
  if (any(vapply(data[vars], function(x) is.numeric(x) && any(!is.finite(x)), logical(1))))
    .stop("Predictors must be finite.")
  # Preserve user functions while making mgcv smooth constructors available.
  env <- new.env(parent = environment(formula))
  for (nm in c("s", "te", "ti", "t2")) env[[nm]] <- get(nm, asNamespace("mgcv"))
  f <- stats::as.formula(call("~", as.name(response), formula[[2]]), env = env)
  f
}
.scale_fit <- function(data, vars) {
  vars <- vars[vapply(data[vars], function(x) is.numeric(x) && length(unique(x)) > 2, logical(1))]
  list(vars = vars, center = vapply(data[vars], mean, numeric(1)),
       scale = vapply(data[vars], sd, numeric(1)))
}
.scale_data <- function(data, scaler) {
  for (v in scaler$vars) data[[v]] <- (data[[v]] - scaler$center[[v]]) / scaler$scale[[v]]
  data
}
.km <- function(data, time, status) {
  fit <- survival::survfit(survival::Surv(data[[time]], 1 - data[[status]]) ~ 1)
  list(time = fit$time, surv = fit$surv)
}
.ipcw <- function(data, time, status, h, km, control, validation = FALSE) {
  t <- data[[time]]; d <- data[[status]]
  j <- findInterval(pmin(t, h), km$time)
  g <- c(1, km$surv)[j + 1L]
  # Historical convention: right-continuous censoring KM, including tied times.
  # Validation weights are not capped (matching the manuscript CV implementation).
  w <- as.numeric(t > h | (t <= h & d == 1)) / pmax(g, control$eps_g)
  if (!validation) w <- pmin(w, control$weight_cap)
  list(y = as.numeric(t <= h & d == 1), w = w, g = g)
}
.both_classes <- function(y, w) {
  if (sum(w * y) <= 0 || sum(w * (1 - y)) <= 0)
    .stop("Each fitting split needs positive IPCW mass for both events and survivors at the horizon.")
}
.with_seed <- function(seed, expr) {
  exists_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (exists_seed) old <- get(".Random.seed", .GlobalEnv)
  on.exit(if (exists_seed) assign(".Random.seed", old, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed); force(expr)
}
.folds <- function(data, time, status, h, folds, seed, fold_id = NULL) {
  n <- nrow(data)
  if (!is.null(fold_id)) {
    if (length(fold_id) != n || anyNA(fold_id)) .stop("fold_id must label every row.")
    f <- match(fold_id, unique(fold_id))
  } else {
    if (length(folds) != 1 || !is.finite(folds) || folds < 2 || folds != as.integer(folds) || folds > n)
      .stop("folds must be an integer between 2 and n.")
    strata <- ifelse(data[[time]] > h, 2L, ifelse(data[[status]] == 1, 1L, 0L))
    f <- .with_seed(seed, {
      a <- integer(n)
      for (s0 in unique(strata)) {
        ix <- which(strata == s0)
        a[ix[sample.int(length(ix))]] <- rep(seq_len(folds), length.out = length(ix))
      }; a
    })
  }
  if (length(unique(f)) < 2 || any(tabulate(f) == 0)) .stop("At least two nonempty folds are required.")
  f
}
.check_gam <- function(fit, label) {
  if (!isTRUE(fit$converged) || any(!is.finite(coef(fit))) || any(!is.finite(fit$sp)))
    .stop(label, " GAM failed to converge; change the specified model or inspect the data.")
  if (!is.null(fit$outer.info$conv) && fit$outer.info$conv != "full convergence")
    .stop(label, " GAM smoothing optimizer: ", fit$outer.info$conv)
}
