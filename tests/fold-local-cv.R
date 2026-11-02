# Exercise the public CV driver and forbid fitting on held-out rows.
library(DTRlearn2)
ns <- asNamespace('DTRlearn2')
cv_data <- getFromNamespace('owl_cv_data', 'DTRlearn2')
fit <- getFromNamespace('owl_fit', 'DTRlearn2')
single <- getFromNamespace('owl_single', 'DTRlearn2')

set.seed(901)
n <- 90L
H <- matrix(rnorm(n*2), ncol=2, dimnames=list(as.character(1:n), NULL))
AA <- lapply(1:3, function(k) sample(c(-1,1), n, replace=TRUE))
RR <- lapply(1:3, function(k) AA[[k]]*(H[,1]+.3*k)+rnorm(n))
known_pi <- rep(list(rep(.5,n)), 3)

.dtr_cv_audit <- new.env()
.dtr_cv_audit$replay <- 0L
.dtr_cv_audit$records <- list()
.dtr_cv_record <- function(kind, x, cost=NULL) {
  stopifnot(.dtr_cv_audit$replay==0L)
  frames <- sys.frames()
  inside <- which(vapply(frames, function(e)
    exists('H_train', envir=e, inherits=FALSE) &&
    exists('H_validation', envir=e, inherits=FALSE), logical(1)))
  in_fold <- length(inside)>0L
  if (in_fold) {
    frame <- frames[[tail(inside,1)]]
    train <- get('H_train', frame)
    validation <- get('H_validation', frame)
    train <- rownames(if (is.list(train)) train[[1]] else train)
    validation <- rownames(if (is.list(validation)) validation[[1]] else validation)
    stopifnot(length(rownames(x))==nrow(x),
      all(rownames(x) %in% train), !any(rownames(x) %in% validation))
  }
  .dtr_cv_audit$records[[length(.dtr_cv_audit$records)+1L]] <-
    list(kind=kind, rows=rownames(x), cost=cost, in_fold=in_fold)
}

run_traced <- function() {
  trace('owl_fit', where=ns, print=FALSE,
    tracer=quote(if (!is.null(cv$fitted)) {
      a <- get('.dtr_cv_audit', .GlobalEnv); a$replay <- a$replay+1L
    }), exit=quote(if (!is.null(cv$fitted)) {
      a <- get('.dtr_cv_audit', .GlobalEnv); a$replay <- a$replay-1L
    }))
  on.exit(untrace('owl_fit', where=ns), add=TRUE)
  trace('owl_residual', where=ns, print=FALSE,
    tracer=quote(get('.dtr_cv_record', .GlobalEnv)('residual', H)))
  on.exit(untrace('owl_residual', where=ns), add=TRUE)
  trace('ql_single', where=ns, print=FALSE,
    tracer=quote(get('.dtr_cv_record', .GlobalEnv)('Q', H)))
  on.exit(untrace('ql_single', where=ns), add=TRUE)
  trace('wsvm_solve', where=ns, print=FALSE,
    tracer=quote(get('.dtr_cv_record', .GlobalEnv)('svm', X, C)))
  on.exit(untrace('wsvm_solve', where=ns), add=TRUE)
  imports <- parent.env(ns)
  trace('cv.glmnet', where=imports, print=FALSE,
    tracer=quote(get('.dtr_cv_record', .GlobalEnv)('glmnet', x)))
  on.exit(untrace('cv.glmnet', where=imports), add=TRUE)

  # Matrix/list histories, one/multiple stages, and both hinge kernels.
  for (K in c(1L,3L)) for (kernel in c('linear','rbf')) {
    for (augment in if (K==1L) FALSE else c(FALSE,TRUE)) {
      set.seed(42)
      model <- owl(H=if (kernel=='linear') H else rep(list(H),K),
        AA=AA[1:K], RR=RR[1:K], n=n, K=K, pi=known_pi[1:K],
        res.lasso=FALSE, kernel=kernel, augment=augment,
        c=c(.5,2), sigma=c(.1,.3), m=3)
      stopifnot(model$type!='has_error', length(model)==K+4L,
        all(vapply(model[1:K], function(x) all(is.finite(x$fit)), logical(1))))
    }
  }
  # Default propensity estimation must also be trained within each fold.
  for (K in 1:2) {
    set.seed(123)
    model <- owl(H=rep(list(H),K), AA=AA[1:K], RR=RR[1:K],
      n=n, K=K, pi='estimated', kernel='rbf', augment=K>1,
      c=c(.5,2), sigma=.1, m=3)
    stopifnot(model$type!='has_error')
  }
  for (kind in c('residual','Q','svm','glmnet'))
    stopifnot(any(vapply(.dtr_cv_audit$records, function(x)
      x$kind==kind && x$in_fold, logical(1))))
  stopifnot(.dtr_cv_audit$replay==0L)
}
run_traced()

# Conditional on selected future parameters, held-out outcomes/features/labels
# cannot affect the training fold's residual weights or learned future models.
set.seed(42)
tuned <- owl(H=rep(list(H),2), AA=AA[1:2], RR=RR[1:2], n=n, K=2,
  pi=known_pi[1:2], augment=TRUE, c=c(.5,2), m=3)
args <- list(H=rep(list(H),2), AA=AA[1:2], RR=RR[1:2], pi=known_pi[1:2],
  n=n, K=2, stage=1, method='hingelinear', pentype='lasso',
  clinear=c(.5,2), sigma=.1, s=1, m=3, solver='svm', tuned=tuned,
  augment=TRUE, folds=rep(1:3, length.out=n))
for (estimate_pi in c(FALSE,TRUE)) {
  args$estimate_pi <- estimate_pi
  set.seed(345); before <- do.call(cv_data, args)
  changed <- args
  held_out <- args$folds==1
  changed$H <- lapply(args$H, function(x) {x[held_out,] <- x[held_out,]+5; x})
  changed$AA <- lapply(args$AA, function(x) {x[held_out] <- -x[held_out]; x})
  changed$RR <- lapply(args$RR, function(x) {x[held_out] <- x[held_out]+100; x})
  set.seed(345); after <- do.call(cv_data, changed)
  stopifnot(identical(before[[1]]$X, after[[1]]$X),
    identical(before[[1]]$Y, after[[1]]$Y), identical(before[[1]]$R, after[[1]]$R),
    !identical(before[[1]]$R_true_t, after[[1]]$R_true_t))
}

# Independent two-stage augmentation formula on held-out observations.
args$estimate_pi <- FALSE
set.seed(345); data <- do.call(cv_data, args)
tr <- args$folds!=1
va <- !tr
center <- colMeans(H[tr,,drop=FALSE])
scale <- apply(H[tr,,drop=FALSE],2,sd)
train_H <- sweep(sweep(H[tr,,drop=FALSE],2,center),2,scale,'/')
x <- sweep(sweep(H[va,,drop=FALSE],2,center),2,scale,'/')
stopifnot(max(abs(data[[1]]$X-train_H))<1e-12,
  max(abs(data[[1]]$Xt-x))<1e-12)
constant_H <- cbind(H[,1], 1)
scaled <- getFromNamespace('owl_scale_fold','DTRlearn2')(
  constant_H[tr,,drop=FALSE], constant_H[va,,drop=FALSE])
stopifnot(all(is.finite(scaled$H)), all(is.finite(scaled$Ht)), all(scaled$H[,2]==0))
set.seed(345)
training <- fit(rep(list(train_H),2), lapply(AA[1:2], function(x) x[tr]),
  lapply(RR[1:2], function(x) x[tr]), sum(tr), 2,
  rep(list(rep(.5,sum(tr))),2), 'hingelinear', 'lasso', TRUE,
  c(.5,2), .1, 1, 3, 'svm', cv=list(target=1, tuned=tuned))
co <- training$qmodels[[2]]$co
Q <- as.numeric(pmax(cbind(1,x,1,x)%*%co, cbind(1,x,-1,-x)%*%co))
future <- training$models[[2]]
d <- getFromNamespace('predict.linearcl','DTRlearn2')(future,x)
expected <- RR[[1]][va]+Q+(RR[[2]][va]-Q)*(as.numeric(d)==AA[[2]][va])/.5
stopifnot(max(abs(data[[1]]$R_true_t-expected))<1e-10)

# Candidate cost remains C/full-supplied-n; validation uses raw augmented rewards.
.dtr_cv_costs <- numeric(0)
check_costs <- function() {
  trace('wsvm_solve', where=ns, print=FALSE, tracer=quote({
    assign('.dtr_cv_costs', c(get('.dtr_cv_costs',.GlobalEnv),C), .GlobalEnv)
  }))
  on.exit(untrace('wsvm_solve',where=ns))
  set.seed(543)
  single(H, AA[[1]], RR[[1]], known_pi[[1]], pentype='LSE',
    clinear=c(.5,2), m=3, cv.data=data)
}
model <- check_costs()
stopifnot(identical(head(.dtr_cv_costs,6), rep(c(.5,2)/n,3)),
  tail(.dtr_cv_costs,1)==model$c/n)
scores <- vapply(c(.5,2), function(C) sum(vapply(data, function(z) {
  candidate <- getFromNamespace('wsvm_solve','DTRlearn2')(
    z$X,z$Y,z$R,C=C/n,solver='svm',e=1e-5)
  d <- getFromNamespace('predict.linearcl','DTRlearn2')(candidate,z$Xt)
  sum(z$R_true_t*(as.numeric(d)==z$Yt)/z$pt)
}, numeric(1))), numeric(1))
stopifnot(model$c==c(.5,2)[which.max(scores)])

# Frozen constant-response Q branches also replay without fitting.
constant <- args
constant$RR[[2]] <- rep(1,n)
set.seed(345); constant_data <- do.call(cv_data, constant)
stopifnot(all(vapply(constant_data, function(z) all(is.finite(z$R_true_t)), logical(1))))
cat('Fold-local CV checks passed, including held-out-row audits and augmentation replay.\n')
