# Compare package scores and treatments with the native WeightSVM model.
library(DTRlearn2)
solve <- getFromNamespace('wsvm_solve', 'DTRlearn2')

check_predictions <- function(model, x, expected) {
  intermediate <- getFromNamespace(paste0('predict.', class(model)), 'DTRlearn2')
  single_input <- if (nrow(x)==1L) as.numeric(x) else x
  stopifnot(identical(as.numeric(intermediate(model, single_input)), expected))
  object <- list(model)
  object$type <- if (inherits(model, 'linearcl')) 'owl_svmlinear' else 'owl_svmrbf'
  class(object) <- 'owl'
  for (full_information in c(FALSE, TRUE)) {
    args <- list(object=object, H=list(x), K=1)
    if (full_information) {
      args$AA <- list(rep(c(-1, 1), length.out=nrow(x)))
      args$RR <- list(rep(1, nrow(x)))
      args$pi <- list(rep(.5, nrow(x)))
    }
    result <- do.call(predict, args)
    stopifnot(identical(as.numeric(result$treatment[[1]]), expected))
  }
}

# Avoid exact zero kernel entries in WeightSVM's precomputed-kernel interface.
X0 <- cbind(c(-3, -2, -1, .5, 2, 4), c(.2, .4, 1.7, -.8, 1.3, -.2))
A0 <- c(-1, -1, 1, -1, 1, 1)
w0 <- c(1, .5, -2, -1, 1.5, .7)
x <- rbind(X0, c(.1, .3), c(1.1, -1.2), c(-2.1, .5))
for (kernel in c('linear', 'rbf')) for (first in c(-1, 1)) {
  first_row <- which(A0 * sign(w0) == first)[1]
  order <- c(first_row, setdiff(seq_len(nrow(X0)), first_row))
  X <- X0[order, , drop=FALSE]
  A <- A0[order]
  wR <- w0[order]
  for (cost in c(.002, 1)) {
    model <- solve(X, A, wR, kernel=kernel, sigma=.3, C=cost, solver='svm')
    gram <- if (kernel=='linear') tcrossprod(X) else kernelMatrix(rbfdot(.3), X)
    cross <- if (kernel=='linear') x %*% t(X) else kernelMatrix(rbfdot(.3), x, X)
    stopifnot(all(gram != 0), all(cross != 0))
    native <- WeightSVM::wsvm(gram, A * sign(wR), weight=abs(wR), cost=cost,
      kernel='precomputed', type='C-classification', scale=FALSE,
      shrinking=FALSE, fitted=FALSE)
    orientation <- as.numeric(native$levels[native$labels[1L]])
    raw <- predict(native, cross, decision.values=TRUE)
    score <- as.numeric(model$beta0 + cross %*% model$alpha1)
    stopifnot(model$svm_orientation == orientation,
      abs(model$beta0 + orientation * as.numeric(native$rho)) < 1e-12,
      max(abs(score - orientation * attr(raw, 'decision.values'))) < 1e-10,
      min(model$alpha1 / (A * sign(wR))) > -1e-10,
      max(model$alpha1 / (A * sign(wR)) - cost * abs(wR)) < 1e-8)
    train_expected <- as.numeric(as.character(predict(native, gram)))
    stopifnot(identical(as.numeric(model$treatment), train_expected))
    expected <- as.numeric(as.character(raw))
    check_predictions(model, x, expected)
    check_predictions(model, x[1, , drop=FALSE], expected[1])
  }
}

# Symmetry gives an exact zero score: native classification picks label two.
for (kernel in c('linear', 'rbf')) for (first in c(-1, 1)) {
  X <- matrix(c(first, -first), ncol=1)
  A <- c(first, -first)
  model <- solve(X, A, c(1, 1), kernel=kernel, sigma=.5, solver='svm')
  origin <- matrix(0, nrow=1)
  cross <- if (kernel=='linear') origin %*% t(X) else kernelMatrix(rbfdot(.5), origin, X)
  # Use raw-feature prediction so the zero margin is independent of that interface.
  native <- WeightSVM::wsvm(X, A, weight=c(1, 1), cost=1,
    kernel=if (kernel=='linear') 'linear' else 'radial', gamma=.5,
    type='C-classification', scale=FALSE,
    shrinking=FALSE, fitted=FALSE)
  stopifnot(as.numeric(model$beta0 + cross %*% model$alpha1) == 0)
  expected <- as.numeric(as.character(predict(native, origin)))
  stopifnot(expected == -model$svm_orientation)
  check_predictions(model, origin, expected)

  # Existing saved models without the new metadata keep their sign() behavior.
  model$svm_orientation <- NULL
  check_predictions(model, origin, 0)
}
cat('Native SVM checks passed: both kernels, class orders, costs and zero scores.\n')
