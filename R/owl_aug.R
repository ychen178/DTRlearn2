
#-----------------------------------------------------------------------------------#
# single-stage augmented outcome-weighted learning (private function)
# Yuan Chen
#-----------------------------------------------------------------------------------#

## Solve single-stage augmented owl
owl_aug <-function(X, AA, RR, n, K, pi, method, pentype='lasso', clinear=2^(-2:2), sigma=c(0.03,0.05,0.07), s=s, m=4, solver='svm',
                  estimate_pi=FALSE, cv=NULL){

  if (K==1) stop(gettextf("No augment methods for single stage data. Please specify augment=F when K=1."))

  select=rep(TRUE,n)
  QL=matrix(0,n,K)
  M=matrix(1,n,K)
  C=matrix(1,n,K)
  models=list()
  qmodels=list()
  prob=matrix(1,n,K)
  QLproj=matrix(0,n,K+1)
  Qspecify=matrix(0,n,K)
  QR_future=0
  Rsum=rep(0, n)
  has_error = 0

  for (k in K:1) {
    A=AA[[k]]

    if (!is.null(cv$fitted)) {
      qmodel = cv$fitted$qmodels[[k]]
      if (is.null(qmodel$co)) {
        QR_future = qmodel$QR_future
        QL[,k] = qmodel$QL
      } else {
        Hk = if (is.list(X)) X[[k]] else X
        QR_future = as.numeric(pmax(cbind(1,Hk,1,Hk) %*% qmodel$co,
          cbind(1,Hk,-1,-Hk) %*% qmodel$co))
        QL[,k] = QR_future
      }
    } else if (min(RR[[k]]+QR_future) != max(RR[[k]]+QR_future)) {
      if(!is.list(X))  output_Q = ql_single(X,A,RR[[k]]+QR_future,pentype=pentype,m=m)
      if(is.list(X))   output_Q = ql_single(X[[k]],A,RR[[k]]+QR_future,pentype=pentype,m=m)
      QR_future = output_Q$Q
      QL[,k] = output_Q$Q
      qmodels[[k]] = list(co=output_Q$co)
    } else {
      QR_future = max(RR[[k]] + QR_future)
      QL[,k] = max(RR[[k]] + QR_future)
      qmodels[[k]] = list(QR_future=QR_future, QL=QL[1,k])
    }
    if(k==K)  R_p = Rsum*select/prob[,K]
    else if (k==K-1)  R_p = Rsum*select/prob[,K] + QLproj[,(k+1):K] * Qspecify[,(k+1):K]
    else  R_p = Rsum*select/prob[,K] + apply(QLproj[,(k+1):K] * Qspecify[,(k+1):K],1,sum)

    R=(RR[[k]]+R_p)
    prob[,k:K] = prob[,k:K] * as.vector(pi[[k]])
    if (method %in% c('hingelinear','hingerbf')) {
      Hk = if (is.list(X)) X[[k]] else X
      if (!is.null(cv) && k==cv$target)
        return(list(H=Hk, A=A, R2=as.numeric(R), pi=pi[[k]], models=models, qmodels=qmodels))
      if (!is.null(cv$fitted)) {
        models[[k]] = cv$fitted$models[[k]]
        predict_stage = if (method=='hingelinear') predict.linearcl else predict.rbfcl
        models[[k]]$treatment = predict_stage(models[[k]], if (nrow(Hk)==1) as.numeric(Hk) else Hk)
      } else {
        cv.data = if (is.null(cv)) owl_cv_data(X, AA, RR, pi, n, K, k, method, pentype,
          clinear, sigma, s, m, solver, models, TRUE, estimate_pi) else NULL
        stage_c = if (is.null(cv)) clinear else cv$tuned[[k]]$c
        stage_sigma = if (is.null(cv)) sigma else cv$tuned[[k]]$sigma
        stage_n = if (is.null(cv)) n else length(cv$tuned[[k]]$fit)
      }
    }

    if(!is.list(X) && is.null(cv$fitted)) {
      if(method == 'hingelinear')
        models[[k]] = owl_single(X,A,R,pi[[k]], pentype=pentype, kernel='linear',clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      if(method == 'hingerbf')
        models[[k]] = owl_single(X,A,R,pi[[k]], pentype=pentype, kernel='rbf',sigma=stage_sigma,clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      if(method %in% c('logit','logitlasso'))
        models[[k]] = owl_logit_single(X,A,R,pi[[k]], pentype=pentype, method=method, m=m)
      if (method %in% c('ols', 'olslasso'))
        models[[k]] = owl_l2_single(X,A,R,pi[[k]], pentype=pentype, method=method, m=m)
      if (method == 'ramplinear')
        models[[k]] = owl_ramp_cv(X,A,R,pi[[k]], pentype=pentype, bigC=clinear, bigS=s, K=m,solver=solver)
    }
    if (is.list(X) && is.null(cv$fitted)) {
      if(method == 'hingelinear')
        models[[k]] = owl_single(X[[k]],A,R,pi[[k]], pentype=pentype, kernel='linear',clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      if(method == 'hingerbf')
        models[[k]] = owl_single(X[[k]],A,R,pi[[k]], pentype=pentype, kernel='rbf',sigma=stage_sigma,clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      if(method %in% c('logit','logitlasso'))
        models[[k]] = owl_logit_single(X[[k]],A,R,pi[[k]], pentype=pentype, method=method, m=m)
      if (method %in% c('ols', 'olslasso'))
        models[[k]] = owl_l2_single(X[[k]],A,R,pi[[k]], pentype=pentype, method=method, m=m)
      if (method == 'ramplinear')
        models[[k]] = owl_ramp_cv(X[[k]],A,R,pi[[k]], pentype=pentype, bigC=clinear, bigS=s, K=m,solver=solver)
    }

    if(is.na(max(models[[k]]$treatment))) {
      has_error = 1
      break
    }

    right = as.vector(models[[k]]$treatment==A)
    select = select * right

    M[,k:K] = M[,k:K] * right
    if (k>1) C[,k:K] = M[,(k-1):(K-1)] - M[,k:K]
    if (k==1) {
      C[,2:K] = M[,1:(K-1)] - M[,2:K]
      C[,1] = rep(1,n) - M[,1]
    }

    Rsum=rep(0, n)
    for (j in k:K){
      if (j>1)  QLproj[,j] = (C[,j]-(1-pi[[j]])*M[,j-1])/prob[,j]
      else      QLproj[,1] = (C[,j]-(1-pi[[j]]))/prob[,j]
      Qspecify[,j] = QL[,j] + Rsum
      Rsum = Rsum + RR[[j]]
    }
  }

  if (method %in% c('hingelinear', 'ramplinear'))  class(models) = 'owl_svmlinear'
  if (method %in% c('hingerbf', 'ramprbf'))  class(models) = 'owl_svmrbf'
  if (method %in% c('logit', 'logitlasso'))  class(models) = 'owl_logit'
  if (method %in% c('ols', 'olslasso'))  class(models) = 'owl_l2'
  if (has_error==1)  class(models) = 'has_error'

  class = class(models)
  res = predict(models, X, AA, RR, K, pi)
  return = c(stage=models, list(valuefun=res$valuefun, benefit=res$benefit, pi=pi, type=class))
  class(return) = "owl"

  return
}





