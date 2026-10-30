library(SPT)
eq <- function(a, b, tol = 1e-8) stopifnot(max(abs(a - b)) < tol)
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
d <- spt_example_data(400, 240, seed = 123)
sf <- ~ sex + s(age, bs = "ts", k = 7) + s(biomarker, bs = "ts", k = 7)
tf <- ~ sex + s(age, bs = "ts", k = 5) + s(biomarker, bs = "ts", k = 5)
src <- spt_source(d$source, "time", "status", sf, "gam", 24)
m <- spt_fit(d$target, "time", "status", 24, tf, src, tau = .5)
stopifnot(m$rho >= 0, m$rho <= 1, m$diagnostics$kkt < 2e-6,
          abs(m$calibration$score) < 1e-10, m$projection_error < 1e-8)
c <- predict(m, d$target, type = "components")
eq(c$target_correction + c$borrowed_projection, c$spt)
eq(predict(m, d$target), pmin(pmax(plogis(c$spt), .001), .999))
f <- tempfile(); saveRDS(m, f); eq(predict(readRDS(f), d$target), predict(m, d$target), 1e-14); unlink(f)

# Check the assembled mgcv penalty against the target anchor score, including
# a two-penalty tensor interaction (not just independent one-penalty smooths).
rich <- ~ sex + s(age, bs="ts", k=5) + s(biomarker, bs="ts", k=5) +
  ti(age, biomarker, bs=c("ts","ts"), k=c(3,3))
ctx <- SPT:::.context(d$target, "time", "status", 24, rich, src, spt_control())
gradient <- crossprod(ctx$B, ctx$w * (plogis(drop(ctx$B %*% coef(ctx$anchor))) - ctx$y)) +
  ctx$S %*% coef(ctx$anchor)
stopifnot(max(abs(gradient)) / sum(ctx$w) < 1e-5)
rich_fit <- spt_fit(d$target, "time", "status", 24, rich, src, tau = 0)
stopifnot(all(is.finite(predict(rich_fit, d$target))))

# Independent joint constrained optimization checks analytic profiling of rho.
ctx <- SPT:::.context(d$target, "time", "status", 24, tf, src, spt_control())
B <- ctx$B; S <- ctx$S; center <- ctx$center; k <- ncol(B)
objective <- function(par) {
  eta <- drop(B %*% par[seq_len(k)]); beta <- par[seq_len(k)] - par[k+1] * center
  b <- pmax(eta,0) + log1p(exp(-abs(eta)))
  (sum((ctx$w+.5)*b - (ctx$w*ctx$y+.5*ctx$q)*eta) + sum(beta*(S%*%beta))/2) / sum(ctx$w)
}
joint <- optim(c(coef(m), m$rho), objective, method="L-BFGS-B",
               lower=c(rep(-Inf,k),0), upper=c(rep(Inf,k),1),
               control=list(maxit=2000, factr=1e3))
eq(joint$value, m$diagnostics$objective, 1e-7)

# CV reconstructs all nuisance fits inside the declared training folds.
ff <- SPT:::.folds(d$target, "time", "status", 24, 3, 13)
cv <- spt_fit(d$target,"time","status",24,tf,src,fold_id=ff)
tr <- d$target[ff!=1,]; va <- d$target[ff==1,]
one <- spt_fit(tr,"time","status",24,tf,src,tau=.5)
a <- SPT:::.ipcw(va,"time","status",24,SPT:::.km(tr,"time","status"),spt_control(),TRUE)
p <- predict(one,va)
eq(cv$cv$losses[ff==1,3],a$w*(-a$y*log(p)-(1-a$y)*log1p(-p)),1e-12)
L <- matrix(c(1,2,1,2, 1.01,1.99,1.01,1.99),4,2)
stopifnot(SPT:::.select_tau(L,c(0,.5),"one_se_half",.5)$tau==.5)
L[,2] <- L[,1]+.2
stopifnot(SPT:::.select_tau(L,c(0,.5),"one_se_half",.5)$tau==0)

# Cox absolute risk agrees with survival's independent survfit interface.
cx <- spt_source(d$source,"time","status",~age+sex+biomarker,"cox")
nd <- d$target[1:5,]
sf0 <- survival::survfit(cx$fit,newdata=SPT:::.scale_data(nd,cx$scaler),se.fit=FALSE)
eq(predict(cx,nd,horizon=24),as.numeric(1-summary(sf0,times=24,extend=TRUE)$surv),1e-12)
eq(plogis(predict(cx,nd,horizon=24,type="link")),
   pmin(pmax(predict(cx,nd,horizon=24),.001),.999),1e-12)
cm <- spt_fit(d$target,"time","status",24,tf,cx,tau=.5)
stopifnot(all(is.finite(predict(cm,nd))))

# Source-only NN tuning, RNG restoration, deterministic fit, and target use.
set.seed(48); before <- .Random.seed
nn <- spt_source(d$source,"time","status",~age+sex+biomarker,"nn",24,
                 folds=2,nn_grid=data.frame(size=c(2L,3L),decay=c(1,1)),nn_starts=2)
stopifnot(identical(before,.Random.seed),sum(nn$tuning$selected)==1)
nm <- spt_fit(d$target,"time","status",24,tf,nn,tau=.5)
stopifnot(all(is.finite(predict(nm,nd))))
f <- tempfile(); saveRDS(nn,f); eq(predict(readRDS(f),nd),predict(nn,nd),1e-14); unlink(f)

# Unpenalized target: rho must be identified by the zero-direction convention.
linear <- spt_fit(d$target,"time","status",24,~sex+age+biomarker,src,tau=.5)
stopifnot(linear$rho==0,linear$diagnostics$zero_penalty_direction)
bad <- d$target; bad$age[1] <- NA
fails(spt_fit(bad,"time","status",24,tf,src,tau=.5))
fails(spt_fit(d$target,"time","status",24,tf,src,tau=-1))
fails(spt_fit(d$target,"time","status",24,tf,src,tau_grid=c(0,1)))
fails(predict(src,nd,horizon=30))
fails(spt_fit(d$target,"time","status",24,~age-1,src,tau=.5))
bad <- d$target; bad$status <- 0
fails(spt_fit(bad,"time","status",24,tf,src,tau=.5))
cat("SPT regression checks passed.\n")
