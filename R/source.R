spt_source <- function(data, time, status, formula, learner = c("gam", "cox", "nn"),
                       horizon = NULL, folds = 5L, fold_id = NULL, seed = 20261029L,
                       nn_grid = expand.grid(size = c(4L, 8L), decay = c(.1, 1, 10)),
                       nn_starts = 3L, control = spt_control()) {
  learner <- match.arg(learner); control <- .control(control)
  .data_check(data, time, status, horizon)
  if (learner != "cox" && is.null(horizon)) .stop("GAM/NN sources require a horizon.")
  f <- .formula(formula, data)
  vars <- all.vars(formula)
  if (any(c(time, status) %in% vars)) .stop("Outcomes cannot appear as predictors.")
  scaler <- .scale_fit(data, vars); z <- .scale_data(data, scaler)
  if (learner == "gam") {
    aux <- .ipcw(data, time, status, horizon, .km(data, time, status), control)
    .both_classes(aux$y, aux$w)
    z$.spt_response <- aux$y; z$.spt_weight <- aux$w
    fit <- mgcv::gam(f, data = z, family = stats::quasibinomial(),
                     weights = .spt_weight, method = "REML", gamma = control$gamma,
                     scale = control$gam_scale, na.action = stats::na.fail)
    .check_gam(fit, "Source")
    obj <- list(fit = fit)
  } else if (learner == "cox") {
    if (!sum(data[[status]])) .stop("Cox source needs events.")
    # Cox uses a conventional linear/factor design, not GAM smooth constructors.
    rhs <- attr(stats::terms(formula), "term.labels")
    if (any(grepl("(^|[^[:alnum:]_])(s|te|ti|t2)\\(", rhs)))
      .stop("Cox source formula must use ordinary model-matrix terms.")
    z$.spt_time <- data[[time]]; z$.spt_event <- data[[status]]
    f[[2]] <- quote(survival::Surv(.spt_time, .spt_event))
    fit <- survival::coxph(f, data = z, ties = "efron", singular.ok = FALSE,
                           model = TRUE, x = TRUE, y = TRUE,
                           control = survival::coxph.control(iter.max = 100, eps = 1e-9),
                           na.action = stats::na.fail)
    if (any(!is.finite(coef(fit))) || any(!is.finite(fit$var)) || max(fit$iter) >= 100)
      .stop("Cox source did not converge to finite coefficients.")
    obj <- list(fit = fit, baseline = survival::basehaz(fit, centered = FALSE),
                max_time = max(data[[time]]))
  } else {
    if (!is.data.frame(nn_grid) || !all(c("size", "decay") %in% names(nn_grid)) ||
        !nrow(nn_grid) || anyNA(nn_grid[c("size", "decay")]) ||
        any(!is.finite(as.matrix(nn_grid[c("size", "decay")])) ) ||
        any(nn_grid$size < 1 | nn_grid$size != as.integer(nn_grid$size)) || any(nn_grid$decay < 0) ||
        length(nn_starts) != 1 || !is.finite(nn_starts) || nn_starts < 1 || nn_starts != as.integer(nn_starts))
      .stop("Invalid NN grid or nn_starts.")
    terms <- stats::terms(formula)
    if (any(grepl("(^|[^[:alnum:]_])(s|te|ti|t2)\\(", attr(terms, "term.labels"))))
      .stop("NN source formula must use ordinary model-matrix terms.")
    ff <- .folds(data, time, status, horizon, folds, seed, fold_id)
    losses <- matrix(NA_real_, nrow(data), nrow(nn_grid))
    for (j in seq_len(nrow(nn_grid))) for (k in sort(unique(ff))) {
      tr <- data[ff != k, , drop = FALSE]; va <- data[ff == k, , drop = FALSE]
      model <- .nn_fit(tr, time, status, horizon, terms, nn_grid[j, ], nn_starts,
                       seed + k * 100L, control)
      p <- .clip(.nn_predict(model, va), control$prob_clip)
      a <- .ipcw(va, time, status, horizon, .km(tr, time, status), control, TRUE)
      losses[ff == k, j] <- a$w * (-a$y * log(p) - (1 - a$y) * log1p(-p))
    }
    best <- which.min(colMeans(losses))
    obj <- .nn_fit(data, time, status, horizon, terms, nn_grid[best, ], nn_starts, seed, control)
    nn_grid$cv_loss <- colMeans(losses); nn_grid$selected <- seq_len(nrow(nn_grid)) == best
    obj$tuning <- nn_grid; obj$fold_id <- ff
    # NN uses fold-specific design-matrix standardization rather than GAM scaling.
    scaler <- NULL
  }
  structure(c(obj, list(learner = learner, horizon = horizon, formula = formula,
                        variables = vars, scaler = scaler, control = control,
                        n_source = nrow(data))), class = "spt_source")
}

.nn_fit <- function(data, time, status, h, terms, grid, starts, seed, control) {
  mf <- stats::model.frame(terms, data, na.action = stats::na.fail)
  mm <- stats::model.matrix(terms, mf)
  cols <- setdiff(colnames(mm), "(Intercept)"); x <- mm[, cols, drop = FALSE]
  center <- colMeans(x); scale <- apply(x, 2, sd)
  if (any(!is.finite(scale) | scale <= 0)) .stop("NN design has constant/invalid columns in a fitting fold.")
  x <- scale(x, center, scale)
  a <- .ipcw(data, time, status, h, .km(data, time, status), control)
  .both_classes(a$y, a$w); w <- a$w * nrow(data) / sum(a$w)
  models <- lapply(seq_len(starts), function(i) .with_seed(seed + i, {
    fit <- nnet::nnet(x = x, y = a$y, weights = w, size = grid$size, entropy = TRUE,
                      linout = FALSE, skip = FALSE, decay = grid$decay,
                      maxit = control$nn_maxit, trace = FALSE, rang = .1,
                      abstol = 1e-8, reltol = 1e-9,
                      MaxNWts = max(2000, (ncol(x) + 2) * grid$size + 1))
    if (fit$convergence != 0 || any(!is.finite(fit$wts))) .stop("NN did not converge.")
    fit
  }))
  list(models = models, terms = terms, columns = cols, center = center, scale = scale,
       contrasts = attr(mm, "contrasts"), xlevels = lapply(mf[vapply(mf, is.factor, logical(1))], levels))
}
.nn_predict <- function(object, newdata) {
  mf <- stats::model.frame(object$terms, newdata, xlev = object$xlevels, na.action = stats::na.fail)
  x <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  x <- scale(x[, object$columns, drop = FALSE], object$center, object$scale)
  rowMeans(vapply(object$models, function(m) as.numeric(predict(m, x, type = "raw")), numeric(nrow(newdata))))
}

predict.spt_source <- function(object, newdata, type = c("response", "link"),
                               horizon = object$horizon, ...) {
  type <- match.arg(type)
  if (is.null(horizon) || length(horizon) != 1 || !is.finite(horizon) || horizon <= 0)
    .stop("Specify one positive horizon.")
  if (object$learner != "cox" && !isTRUE(all.equal(horizon, object$horizon)))
    .stop("GAM/NN sources are horizon-specific; fit another source at this horizon.")
  .formula(object$formula, newdata)
  if (object$learner == "gam") {
    eta <- as.numeric(predict(object$fit, .scale_data(newdata, object$scaler), type = "link"))
    if (any(!is.finite(eta))) .stop("Non-finite source GAM prediction.")
    # Preserve the historical raw GAM logit score for projection/calibration.
    return(if (type == "link") eta else plogis(eta))
  }
  if (object$learner == "nn") p <- .nn_predict(object, newdata) else {
    if (horizon > object$max_time) .stop("Cox horizon exceeds source observed follow-up.")
    lp <- as.numeric(predict(object$fit, .scale_data(newdata, object$scaler),
                              type = "lp", reference = "zero"))
    j <- findInterval(horizon, object$baseline$time)
    H <- if (j == 0L) 0 else object$baseline$hazard[j]
    p <- if (H == 0) rep(0, length(lp)) else -expm1(-exp(log(H) + lp))
  }
  if (any(!is.finite(p))) .stop("Non-finite source prediction.")
  if (type == "link") qlogis(.clip(p, object$control$prob_clip)) else p
}
print.spt_source <- function(x, ...) {
  cat("SPT source:", x$learner, "| n =", x$n_source, "| horizon =",
      if (is.null(x$horizon)) "supplied at prediction" else x$horizon, "\n")
  invisible(x)
}
