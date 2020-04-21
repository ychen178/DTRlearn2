
#-----------------------------------------------------------------------------------#
# single-stage (private function) and general K-stage Q-learning
# Yuan Chen, April 2020
#-----------------------------------------------------------------------------------#

## solve a single-stage Q-learning
ql_single <- function(H, A, R, pentype='lasso', m=4){
  n = length(A)
  X = cbind(H, A, diag(A)%*%H)
  cvfit = cv.glmnet(X, R, nfolds=m)
  
  if (pentype=='lasso'){
    co = as.matrix(predict(cvfit, s="lambda.min", type="coeff"))
  } else {
    co = predict(cvfit, s=0, type="coeff")
  }
  
  XX1 = cbind(rep(1,n), H, rep(1,n), diag(n)%*%H)
  XX2 = cbind(rep(1,n), H, rep(-1,n),-diag(n)%*%H)
  Q1 = XX1 %*% co
  Q2 = XX2 %*% co
  Q = apply(cbind(Q1, Q2), 1, max)
  
  # no tailoring variables chosen (interaction terms), randomize the recommneded treatment
  treatment = 2 * (Q1 > Q2) - 1
  if (sum(Q1==Q2)==n)  treatment = rbinom(n, 1, 0.5)
  
  Qsingle = list(co=co, Q=Q, treatment = treatment)
  class(Qsingle) = 'qsingle'
  Qsingle
}


### a general K-stage Q-learning
ql <- function (H, AA, RR, K, pi='estimated', lasso=TRUE, m=4) {
  
  if(lasso==TRUE) pentype = 'lasso'
  else if (lasso==FALSE) pentype = 'LSE'
  else print("lasso should be TRUE or FALSE")
  
  if(pi[[1]][1]=='estimated' & K>1) {
    pi = list()
    for (j in 1:K) {
      Y = AA[[j]]*0.5 + 0.5
      if (is.list(H)) Hj = H[[j]]
      else Hj = H
      pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=m)
      pred = predict(pi_logit, Hj, s="lambda.min")
      pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
    }
  }
  if(pi[[1]][1]=='estimated' & K==1) {
    Y = AA*0.5 + 0.5
    pi_logit = cv.glmnet(H, Y, family = "binomial", nfolds=m)
    pred = predict(pi_logit, H, s="lambda.min")
    pi = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
  }
  
  if (K==1) {
    if (!is.list(AA)) AA = list(AA)
    if (!is.list(RR)) RR = list(RR)
    if (!is.list(pi)) pi = list(pi)
  }
  
  R_future=0
  coef = list()
  results = list()
  
  if (! is.list(H)) {
    for (j in K:1) {
      R=RR[[j]]+R_future
      if (min(R)!=max(R)) {
        results[[j]] = ql_single(H, AA[[j]], R, pentype=pentype, m=m)
        R_future = results[[j]]$Q
      } else {
        results[[j]]=list(co=rep(0,2+2*dim(H)[2]),Q=R)
        R_future=R
      }
    }
  }
  if (is.list(H)) {
    for (j in K:1) {
      R=RR[[j]]+R_future
      if (min(R) != max(R)) {
        results[[j]] = ql_single (H[[j]], AA[[j]], R, pentype=pentype, m=m)
        R_future = results[[j]]$Q
      } else {
        results[[j]] = list(co = rep(0,2+2*dim(H[[j]])[2]), Q = R)
        R_future = R
      }
    }
  }
  class(results) = 'ql'
  res = predict(results, H, AA, RR, K, pi)
  return = c(stage=results, list(valuefun=res$valuefun, benefit=res$benefit, pi=pi))
  class(return) = 'ql'
  return
}





