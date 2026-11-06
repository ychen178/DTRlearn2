# Private helpers for fold-local hinge-loss tuning.
owl_residual <- function(H, R2, pi, pentype, m) {
  n = length(R2)
  if (max(R2)!=min(R2)) {
    if (pentype=='lasso') {
      cvfit = cv.glmnet(H, R2, nfolds=m)
      co = as.matrix(predict(cvfit, s='lambda.min', type='coeff'))
    }
    if (pentype=='LSE') co = coef(lm(R2~H))
    r = R2 - cbind(rep(1,n), H) %*% co
  } else r = R2
  r / pi
}

owl_folds <- function(n, m) {
  rand = sample(m, n, replace=TRUE)
  for (attempt in 1:10) {
    counts = sort(tabulate(rand, nbins=m))
    if (counts[1]>=2 && !any(counts[1:(m-1)]<3)) return(rand)
    if (attempt<10) rand = sample(m, n, replace=TRUE)
  }
  stop('Too few observations in the training or test set for cross validation. Please increase sample size.')
}

# H is standardized by the caller; restandardize CV folds using training rows.
owl_scale_fold <- function(H, Ht) {
  if (is.list(H)) {
    scaled = Map(owl_scale_fold, H, Ht)
    return(list(H=lapply(scaled, function(x) x$H), Ht=lapply(scaled, function(x) x$Ht)))
  }
  center = colMeans(H)
  scale = apply(H, 2, stats::sd)
  scale[!is.finite(scale) | scale==0] = 1
  list(H=sweep(sweep(H,2,center),2,scale,'/'),
    Ht=sweep(sweep(Ht,2,center),2,scale,'/'))
}

# Refit the already selected future stages on training rows only. Replay the
# original backward recursion on validation rows with those fitted models.
owl_cv_data <- function(H, AA, RR, pi, n, K, stage, method, pentype,
                        clinear, sigma, s, m, solver, tuned, augment,
                        estimate_pi=FALSE, folds=owl_folds(n,m)) {
  subset_H <- function(rows) {
    if (is.list(H)) lapply(H, function(x) x[rows,,drop=FALSE])
    else H[rows,,drop=FALSE]
  }
  data = vector('list', m)
  for (i in 1:m) {
    train = folds!=i
    validation = !train
    H_train = subset_H(train)
    H_validation = subset_H(validation)
    scaled = owl_scale_fold(H_train, H_validation)
    H_train = scaled$H
    H_validation = scaled$Ht
    AA_train = lapply(AA, function(x) x[train])
    AA_validation = lapply(AA, function(x) x[validation])
    pi_train = lapply(pi, function(x) x[train])
    pi_validation = lapply(pi, function(x) x[validation])
    if (estimate_pi) for (j in 1:K) {
      X = if (is.list(H_train)) H_train[[j]] else H_train
      Xt = if (is.list(H_validation)) H_validation[[j]] else H_validation
      model = cv.glmnet(X, (AA_train[[j]]+1)/2, family='binomial', nfolds=3)
      p = as.numeric(predict(model, X, s='lambda.min', type='response'))
      pt = as.numeric(predict(model, Xt, s='lambda.min', type='response'))
      pi_train[[j]] = ifelse(AA_train[[j]]==1, p, 1-p)
      pi_validation[[j]] = ifelse(AA_validation[[j]]==1, pt, 1-pt)
    }
    training = owl_fit(H_train, AA_train, lapply(RR, function(x) x[train]),
      sum(train), K, pi_train, method, pentype, augment, clinear, sigma, s,
      m, solver, cv=list(target=stage, tuned=tuned))
    if (is.null(training$H)) stop('Too few training observations remain after future-stage selection in cross validation.')
    testing = owl_fit(H_validation, AA_validation, lapply(RR, function(x) x[validation]),
      sum(validation), K, pi_validation, method, pentype, augment, clinear, sigma, s,
      m, solver, cv=list(target=stage, fitted=training))
    data[[i]] = list(X=training$H, Y=training$A,
      R=owl_residual(training$H, training$R2, training$pi, pentype, m),
      Xt=testing$H, Yt=testing$A, R_true_t=testing$R2, pt=testing$pi)
  }
  data
}
