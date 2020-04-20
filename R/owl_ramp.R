
#-----------------------------------------------------------------------------------#
# functions for outcome-weighted learning with ramp loss (private functions)
# Yuan Chen, April 2020
#-----------------------------------------------------------------------------------#

ramp.loss <- function(u,s){
  (u/s <= -1/2) + (u/s > -1/2) * (u/s < 1/2) * (-u/s + 1/2)
}


### solve single-stage owl with ramp loss under a given set of c and s
owl_ramp_notune = function (X, A, R, pi, newX, newA, newR, newpi, pentype='lasso', iter.Max=20, beta.init=NULL, tol=1e-5, C=0.1, s=0.1) {

  n = length(A)
  if (max(R) != min(R)) {
    if (pentype=='lasso') {
      cvfit = cv.glmnet(X, R)
      co = as.matrix(predict(cvfit, s="lambda.min", type="coeff"))
    }
    if (pentype == 'LSE'){
      co = coef(lm(R ~ X))
    }
    W = R - cbind(rep(1,n), X) %*% co
  }
  else W = R

  Wabs <- abs(W)/pi
  Z <- (W>0) - (W<=0)
  As <- Z*A

  Q <- X%*%t(X) * (As%*%t(As))
  Q1 <- X%*%t(X) * ((As*Wabs)%*%t(As*Wabs))
  p <- dim(X)[2]
  delta <- 99999

  # initial values of beta from owl under hinge loss
  if(is.null(beta.init)){
    fit_temp <- owl_single(X, A, R, pi, kernel="linear")
    beta.init <- c(fit_temp$beta0, fit_temp$beta)
  }

  eta.init <- As * cbind(rep(1, n), X) %*% beta.init
  ramploss.init <- C * sum(Wabs*ramp.loss(eta.init, s)) + 0.5 * sum (beta.init[-1]^2)
  stop_loop <- 0
  iter <- 1
  has_error = 0

  #iteratively solve qudartive programming unitl ramp.loss converges
  while (iter < iter.Max ){
    dh2 = -(eta.init < -0.5 * s)/s
    f = -(1 - 2*C*Q %*% (Wabs* dh2)/s)/2 *Wabs

    H = Q1/s^2
    H = H + 1e-8 * diag(NCOL(K)) %*% (tcrossprod(Wabs))
    fit <- tryCatch(ipop(f, H, A=t(As*Wabs)/s, b=-C*sum(As*Wabs*dh2), l=rep(0, n), u=rep(C,n), r=0), error=function(e) e)

    if("error" %in% class(fit))  {
      has_error = 1
      print(fit)
      break
    }

    alpha1 <- primal(fit)
    alpha <- alpha1 * Wabs
    beta.pred <- t(X) %*% (As* (C*Wabs*dh2 + alpha/s))

    beta0 <- s/As/2 - X%*% beta.pred
    index <- (alpha1 >C * 1e-4) & (alpha1 < C*(1-1e-4))
    if (sum(index) > 0){
      beta0.pred <- mean(beta0[index])
    } else {
      beta0.pred <- -max((As==-1)*beta0)/2 - min((As==1)*beta0)/2
    }

    beta.new <- c(beta0.pred, beta.pred)
    eta.new <- As * cbind(rep(1, n), X)%*% beta.new
    ramploss.new <-  C * sum(Wabs*ramp.loss(eta.new, s)) + 0.5 * sum (beta.new[-1]^2)

    delta <- abs(ramploss.new-ramploss.init)

    if (delta < tol) break

    iter <- iter + 1
    beta.init <- beta.new
    eta.init <- eta.new
    ramploss.init <- ramploss.new
  }

  if (has_error == 0) {
    fit = cbind(rep(1, nrow(newX)), newX)%*% beta.new
    prob = exp(fit) / (1+ exp(fit))
    treatment = (fit>0)*2-1
    valuefun = mean(newR*(treatment==newA)/newpi)
    benefit = valuefun - mean(newR*(treatment!=newA)/newpi)
  }
  else {
    beta0.pred = NA
    beta.pred = NA
    fit = NA
    prob=NA
    treatment = NA
    valuefun = NA
    benefit = NA
  }
  return(list(beta0=beta0.pred, beta=beta.pred, fit=fit, probability=prob, treatment=treatment, valuefun = valuefun, benefit = benefit, iter=iter))
}


## Solve single-stage owl with ramp loss and CV to choose tuning parameters c and s
owl_ramp_cv  <- function(X, A, Y, PrTx, pentype=pentype, bigC, bigS, K, iter.Max=20, tol=1e-5, beta.init=NULL) {

  n = length(A)
  bigC_input = bigC
  bigC = bigC / n

  folds = sample(1:K, n, replace=T)
  folds_n = 1
  while (folds_n < 10) {
    folds_nobs = rep(0, K)
    for (t in 1:K)  folds_nobs[t] = sum(folds==t)
    folds_nobs = sort(folds_nobs)

    if(folds_nobs[1]< 2 | sum(folds_nobs[1:(K-1)]< 3))   folds = sample(1:K, n, replace=T)
    else  break
    folds_n = folds_n + 1
  }
  if (folds_n == 10) {
    print("Too few observations in the training or test set for cross validation. Please increase sample size.")
    return(list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, valuefun=NA, benefit=NA, c=NA, s=NA, iter=NA))
  }

  Vtemp <- -99999
  OC = NA
  OS = NA

  for(i in 1:length(bigC)){
    if(is.null(beta.init)){
      fit_temp <- owl_single(X, A, Y, PrTx, kernel = "linear")
      beta.init = c(fit_temp$beta0, fit_temp$beta)
      if (is.na(max(beta.init)))  {
        beta.init = rnorm(dim(X)[2]+1, 0, 1)
      }
    }
    for (j in 1:length(bigS)) {
      predV <- rep(NA, K)
      for (t in 1:K) {
        Ytrain <- Y[folds!=t]; Xtrain <- X[folds!=t,]; Atrain <- A[folds!=t]; PTtrain <- PrTx[folds!=t];
        Ytest <- Y[folds==t]; Xtest <- X[folds==t,];  Atest <- A[folds==t]; PTtest <- PrTx[folds==t];
        temp = owl_ramp_notune(Xtrain, Atrain, Ytrain, PTtrain, Xtest, Atest, Ytest, PTtest, pentype=pentype, C=bigC[i], s=bigS[j], beta.init = beta.init, tol=tol)
        predV[t] <- sum(Ytest*(Atest==temp$treatment)/PTtest)
      }
      if (!is.na(sum(predV)) & sum(predV)>Vtemp*length(A)) {
        OC=i; OS=j
        Vtemp <- sum(predV)/length(A)
      }
    }
  }
  if( !is.na(OC) & !is.na(OS) ) {
    model = owl_ramp_notune(X, A, Y, PrTx, X, A, Y, PrTx, pentype=pentype, iter.Max=iter.Max, beta.init=beta.init, tol=tol, C=bigC[OC], s=bigS[OS])
    result = c(model, c=bigC_input[OC], s=bigS[OS])
  }
  else result = list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, valuefun=NA, benefit=NA, c=NA, s=NA, iter=NA)
  return(result)
}










