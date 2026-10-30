.penalty <- function(fit) {
  if (!is.null(fit$full.sp)) .stop("Linked smoothing parameters (smooth id) are not supported in this release.")
  S <- matrix(0, length(coef(fit)), length(coef(fit)))
  for (sm in fit$smooth) {
    ix <- sm$first.para:sm$last.para
    if (length(sm$S)) for (j in seq_along(sm$S)) {
      sp <- fit$sp[sm$first.sp + j - 1L]
      if (length(sp) != 1 || !is.finite(sp) || sp < 0) .stop("Invalid smooth penalty mapping.")
      S[ix, ix] <- S[ix, ix] + sp * sm$S[[j]]
    }
  }
  (S + t(S)) / 2
}
.context <- function(data, time, status, h, formula, source, control) {
  vars <- all.vars(formula); scaler <- .scale_fit(data, vars)
  z <- .scale_data(data, scaler); f <- .formula(formula, z)
  km <- .km(data, time, status); a <- .ipcw(data, time, status, h, km, control)
  .both_classes(a$y, a$w)
  z$.spt_response <- a$y; z$.spt_weight <- a$w
  anchor <- mgcv::gam(f, data = z, weights = .spt_weight, family = stats::quasibinomial(),
                      method = "REML", gamma = control$gamma, scale = control$gam_scale,
                      na.action = stats::na.fail)
  .check_gam(anchor, "Target")
  if (!is.null(anchor$offset) && any(anchor$offset != 0)) .stop("Target offsets are not supported.")
  B <- unname(predict(anchor, z, type = "lpmatrix")); qb <- qr(B, tol = 1e-10)
  if (qb$rank != ncol(B)) .stop("Target basis is not full column rank; simplify the user-specified formula.")
  s <- predict(source, data, type = "link", horizon = h)
  z$.spt_response <- s
  projection <- mgcv::gam(f, data = z, family = stats::gaussian(), method = "REML",
                          gamma = control$gamma, na.action = stats::na.fail)
  .check_gam(projection, "Projection")
  projected <- as.numeric(predict(projection, z, type = "response"))
  center <- drop(qr.coef(qb, projected))
  error <- max(abs(drop(B %*% center) - projected))
  if (error > 1e-8) .stop("Source projection cannot be represented in target basis.")
  score <- function(alpha) sum(a$w * (plogis(s + alpha) - a$y)) / sum(a$w)
  lo <- -8; hi <- 8
  while (score(lo) > 0 && lo > -1024) lo <- lo * 2
  while (score(hi) < 0 && hi < 1024) hi <- hi * 2
  alpha <- stats::uniroot(score, c(lo, hi), tol = 1e-12)$root
  q <- .clip(plogis(s + alpha), control$prob_clip)
  list(anchor = anchor, projection = projection, scaler = scaler, B = B,
       center = center, S = .penalty(anchor), y = a$y, w = a$w, q = q, km = km,
       calibration = list(alpha = alpha, score = score(alpha),
                          clipped_score = sum(a$w * (q - a$y)) / sum(a$w)),
       projection_error = error, residual_rms = sqrt(mean((s - projected)^2)))
}

.solve <- function(B, center, y, w, q, tau, S, start, tolerance) {
  k <- ncol(B); v <- w + tau; wy <- w * y + tau * q; norm <- sum(w)
  Sc <- drop(S %*% center); den <- sum(center * Sc)
  zero <- abs(den) <= 1e-14 * max(1, sum(abs(center * Sc)))
  if (den < -1e-8) .stop("Negative source penalty curvature.")
  R <- chol(crossprod(B, B * (v / 4)) + S); Q <- backsolve(R, diag(k))
  unpack <- function(z) {
    theta <- drop(Q %*% z)
    rho <- if (zero) 0 else min(1, max(0, sum(theta * Sc) / den))
    list(theta = theta, rho = rho, beta = theta - rho * center, eta = drop(B %*% theta))
  }
  value <- function(z) {
    a <- unpack(z)
    (sum(v * .softplus(a$eta) - wy * a$eta) + sum(a$beta * (S %*% a$beta)) / 2) / norm
  }
  gradient <- function(z) {
    a <- unpack(z)
    drop(crossprod(Q, crossprod(B, v * plogis(a$eta) - wy) + S %*% a$beta)) / norm
  }
  opt <- stats::optim(drop(R %*% start), value, gradient, method = "L-BFGS-B",
                      control = list(maxit = 3000, factr = 1e3, pgtol = 1e-10, lmm = 20))
  attempts <- list(list(par = opt$par, value = value(opt$par), code = opt$convergence,
                        kkt = max(abs(gradient(opt$par))), method = "L-BFGS-B"))
  if (opt$convergence != 0 || attempts[[1]]$kkt > tolerance) {
    opt <- stats::nlminb(opt$par, value, gradient,
                         control = list(iter.max = 10000, eval.max = 20000, rel.tol = 1e-12, x.tol = 1e-12))
    attempts[[2]] <- list(par = opt$par, value = value(opt$par), code = opt$convergence,
                          kkt = max(abs(gradient(opt$par))), method = "nlminb")
  }
  ok <- vapply(attempts, function(a) a$code == 0 && a$kkt <= tolerance, logical(1))
  values <- vapply(attempts, `[[`, numeric(1), "value")
  if (!any(ok)) .stop("SPT optimizer failed convergence/stationarity checks.")
  j <- which(ok)[which.min(values[ok])]; best <- attempts[[j]]
  if (best$value - min(values) > 1e-10) .stop("Converged SPT solution has an inferior objective.")
  if (length(attempts) > 1) warning("SPT used an unchanged-objective optimizer restart.", call. = FALSE)
  a <- unpack(best$par); g <- -sum(center * (S %*% a$beta)) / norm
  rho_kkt <- abs(if (a$rho <= 1e-8) min(g, 0) else if (a$rho >= 1 - 1e-8) max(g, 0) else g)
  score <- sum(v * plogis(a$eta) - wy) / (sum(w) + tau * nrow(B))
  if (rho_kkt > tolerance || abs(score) > 1e-5) .stop("SPT profile/intercept score check failed.")
  list(theta = a$theta, beta = a$beta, rho = a$rho, tau = tau,
       diagnostics = list(objective = best$value, kkt = best$kkt, rho_kkt = rho_kkt,
                          intercept_score = score, zero_penalty_direction = zero,
                          rho_penalty_curvature = den / norm, optimizer = best$method),
       attempts = lapply(attempts, function(a) a[setdiff(names(a), "par")]))
}
.fit_tau <- function(ctx, tau, control) {
  .solve(ctx$B, ctx$center, ctx$y, ctx$w, ctx$q, tau, ctx$S,
         coef(ctx$anchor), control$kkt_tol)
}
.select_tau <- function(losses, grid, rule, reference) {
  risk <- colMeans(losses); best <- order(risk, grid)[1]
  differences <- sweep(losses, 1, losses[, best], "-")
  delta <- colMeans(differences); se <- apply(differences, 2, sd) / sqrt(nrow(losses))
  eligible <- delta <= se + 64 * .Machine$double.eps * max(1, max(abs(risk)))
  chosen <- best
  if (rule == "one_se_half") {
    ref <- which(grid == reference)
    if (eligible[ref]) chosen <- ref
  }
  list(tau = grid[chosen], table = data.frame(tau = grid, cv_loss = risk,
       difference_to_best = delta, paired_se = se, eligible = eligible,
       empirical_best = seq_along(grid) == best, selected = seq_along(grid) == chosen))
}

spt_fit <- function(data, time, status, horizon, target_formula, source,
                    tau = NULL, tau_grid = c(0, .1, .5, 1),
                    selection = c("one_se_half", "min_cv"), reference_tau = .5,
                    folds = 3L, fold_id = NULL, seed = 20261029L,
                    control = spt_control()) {
  control <- .control(control); selection <- match.arg(selection)
  .data_check(data, time, status, horizon); .formula(target_formula, data)
  if (!inherits(source, "spt_source")) .stop("source must be fitted by spt_source().")
  if (any(c(time, status) %in% all.vars(target_formula))) .stop("Outcomes cannot be target predictors.")
  # Require an intercept so calibration and the intercept score have their stated meaning.
  if (attr(stats::terms(target_formula), "intercept") != 1L ||
      length(attr(stats::terms(target_formula), "offset")))
    .stop("Target formula must include an intercept and cannot include offsets.")
  cv <- NULL; ff <- NULL
  if (is.null(tau)) {
    if (!is.numeric(tau_grid) || !length(tau_grid) || any(!is.finite(tau_grid)) ||
        any(tau_grid < 0) || anyDuplicated(tau_grid)) .stop("tau_grid must contain unique finite nonnegative values.")
    if (selection == "one_se_half" && (length(reference_tau) != 1 ||
        !is.finite(reference_tau) || sum(tau_grid == reference_tau) != 1))
      .stop("reference_tau must appear exactly once in tau_grid.")
    ff <- .folds(data, time, status, horizon, folds, seed, fold_id)
    losses <- matrix(NA_real_, nrow(data), length(tau_grid))
    for (f in sort(unique(ff))) {
      tr <- data[ff != f, , drop = FALSE]; va <- data[ff == f, , drop = FALSE]
      ctx <- .context(tr, time, status, horizon, target_formula, source, control)
      V <- predict(ctx$anchor, .scale_data(va, ctx$scaler), type = "lpmatrix")
      a <- .ipcw(va, time, status, horizon, ctx$km, control, TRUE)
      for (j in seq_along(tau_grid)) {
        model <- .fit_tau(ctx, tau_grid[j], control)
        p <- .clip(plogis(drop(V %*% model$theta)), control$prob_clip)
        losses[ff == f, j] <- a$w * (-a$y * log(p) - (1 - a$y) * log1p(-p))
      }
    }
    cv <- .select_tau(losses, tau_grid, selection, reference_tau)
    tau <- cv$tau; cv$losses <- losses; cv$fold_id <- ff
  } else if (length(tau) != 1 || !is.finite(tau) || tau < 0) .stop("tau must be NULL or one finite nonnegative number.")
  ctx <- .context(data, time, status, horizon, target_formula, source, control)
  fit <- .fit_tau(ctx, tau, control)
  names(fit$theta) <- names(fit$beta) <- names(ctx$center) <- names(coef(ctx$anchor))
  structure(c(fit, list(anchor = ctx$anchor, projection = ctx$projection,
                        projection_coef = ctx$center, penalty = ctx$S, scaler = ctx$scaler,
                        calibration = ctx$calibration, projection_error = ctx$projection_error,
                        residual_rms = ctx$residual_rms, horizon = horizon,
                        formula = target_formula, n_target = nrow(data),
                        source_learner = source$learner, cv = cv,
                        selection = if (is.null(cv)) "fixed" else selection,
                        control = control, call = match.call())), class = "spt")
}

predict.spt <- function(object, newdata, type = c("response", "link", "components"), ...) {
  type <- match.arg(type); .formula(object$formula, newdata)
  B <- predict(object$anchor, .scale_data(newdata, object$scaler), type = "lpmatrix")
  eta <- drop(B %*% object$theta)
  if (any(!is.finite(eta))) .stop("Non-finite target predictions; inspect newdata and factor levels.")
  if (type == "link") return(eta)
  if (type == "response") return(.clip(plogis(eta), object$control$prob_clip))
  data.frame(target = drop(B %*% coef(object$anchor)),
             projection = drop(B %*% object$projection_coef),
             target_correction = drop(B %*% object$beta),
             borrowed_projection = object$rho * drop(B %*% object$projection_coef),
             spt = eta)
}
coef.spt <- function(object, ...) object$theta
print.spt <- function(x, ...) {
  cat("Source Projection Transfer (SPT)\n", "Source:", x$source_learner,
      "| target n:", x$n_target, "| horizon:", x$horizon, "\n",
      "tau:", x$tau, "(", x$selection, ") | rho:", format(x$rho, digits = 4), "\n",
      "Calibration intercept:", format(x$calibration$alpha, digits = 4),
      "| KKT:", format(x$diagnostics$kkt, digits = 3), "\n")
  invisible(x)
}
