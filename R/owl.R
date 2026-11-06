
#-----------------------------------------------------------------------------------#
# functions for the single-stage (private) and general K-stage outcome-weighted learning
# Yuan Chen
#-----------------------------------------------------------------------------------#

## solve single stage outcome-weighted learning
owl_single <-function(H, A, R2, pi, pentype='lasso', kernel='linear', sigma=c(0.03,0.05,0.07), clinear=2.^(-2:2), m=3, e=1e-5, solver='svm',
                      cv.data=NULL, tune=TRUE, n_total=length(A))  {
  
  npar=length(clinear)
  n=length(A)
  p=dim(H)[2]
  clinear_input = clinear
  clinear = clinear/n_total
  
  r = owl_residual(H, R2, pi, pentype, m)
  if (!tune) {
    model = wsvm_solve(H, A, r, kernel=kernel, sigma=sigma[1], C=clinear[1], e=e, solver=solver)
    if (kernel=='linear') return(c(model, c=clinear_input[1]))
    return(c(model, c=clinear_input[1], sigma=sigma[1]))
  }

  if (is.null(cv.data)) {
    rand = owl_folds(n,m)
    cv.data = lapply(1:m, function(i) {
      this = rand!=i
      scaled = owl_scale_fold(H[this,,drop=FALSE], H[!this,,drop=FALSE])
      list(X=scaled$H, Y=A[this],
        R=owl_residual(scaled$H, R2[this], pi[this], pentype, m),
        Xt=scaled$Ht, Yt=A[!this], R_true_t=R2[!this], pt=pi[!this])
    })
  }
  
  if (kernel=='linear'){
    V=matrix(NA,m,npar)
    for (i in 1:m) {
      X=cv.data[[i]]$X
      Y=cv.data[[i]]$Y
      R=cv.data[[i]]$R
      Xt=cv.data[[i]]$Xt
      Yt=cv.data[[i]]$Yt
      R_true_t=cv.data[[i]]$R_true_t
      pt=cv.data[[i]]$pt
      for (j in 1:npar) {
        model = wsvm_solve(X, Y, R, C=clinear[j], e=e, solver=solver)
        if (! is.na(max(model$alpha1)) ) {
          YP = if (nrow(Xt)) predict(model, Xt) else numeric(0)
          V[i,j] = sum(R_true_t*(YP==Yt)/pt)
        }
      }
    }
    mimi = apply(V, 2, sum)
    if (max(mimi, na.rm=T) == -Inf)
      result = list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, c=NA, alpha1=NA)
    else {
      best=which.max(mimi)
      cbest=clinear[best]
      model=wsvm_solve(H,A,r,C=cbest,e=e,solver=solver)
      result = c(model, c=clinear_input[best])
    }
  }
  
  if (kernel=='rbf') {
    nsig=length(sigma)
    V=array(0,c(npar,nsig,m))
    for (i in 1:m){
      X=cv.data[[i]]$X
      Y=cv.data[[i]]$Y
      R=cv.data[[i]]$R
      Xt=cv.data[[i]]$Xt
      Yt=cv.data[[i]]$Yt
      R_true_t=cv.data[[i]]$R_true_t
      pt=cv.data[[i]]$pt
      for (j in 1:npar){
        for (s in 1:nsig){
          model=wsvm_solve(X,Y,R,'rbf',sigma=sigma[s],C=clinear[j],e=e,solver=solver)
          if (! is.na(max(model$alpha1)) ) {
            YP=if (nrow(Xt)) predict(model, if (nrow(Xt)==1) as.numeric(Xt) else Xt) else numeric(0)
            V[j,s,i]=sum(R_true_t*(YP==Yt)/pt)
          }
        }
      }
    }
    mimi=apply(V,c(1,2),sum)
    if (max(mimi, na.rm=T) == -Inf)
      result = list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, c=NA, sigma=NA, H=NA, alpha1=NA)
    else {
      best=which(mimi==max(mimi, na.rm=T),arr.ind=TRUE)
      bestC=clinear[best[1,1]]
      bestSig=sigma[best[1,2]]
      model = wsvm_solve(H,A,r,'rbf', bestSig, C=bestC, e=e, solver=solver)
      result = c(model, c=clinear_input[best[1,1]], sigma=bestSig)
    }
  }
  result
}



## a general K-stage outcome-weighted learning (changed default to augment=FALSE)
owl <- function(H,AA,RR,n,K,pi='estimated', res.lasso=TRUE, loss='hinge', kernel='linear', augment=FALSE, c=2.^(-2:2), sigma=c(0.03,0.05,0.07), s=2.^(-2:2), m=4, solver='svm') {
  
  if(res.lasso==TRUE) pentype = 'lasso'
  else if (res.lasso==F) pentype = 'LSE'
  else print("res.lasso should be TRUE or FALSE")
  
  if(loss=='hinge') {
    if (kernel=='linear') method = 'hingelinear'
    else if (kernel=='rbf') method = 'hingerbf'
    else stop(gettextf("Please specify kernel to be 'linear' or 'rbf' for hinge loss."))
  }
  else if (loss=='ramp') {
    if (kernel=='linear') method = 'ramplinear'
    else stop(gettextf("Please specify kernel to be 'linear' for ramp loss."))
  }
  else if (loss=='logit.lasso') method = 'logitlasso'
  else if (loss=='logit') method = 'logit'
  else if (loss=='l2.lasso') method = 'olslasso'
  else if (loss=='l2') method = 'ols'
  else stop(gettextf("Please specify a valid loss function."))
  
  estimate_pi = pi[[1]][1]=='estimated'
  if(estimate_pi & K > 1) {
    pi = list()
    for (j in 1:K) {
      Y = AA[[j]]*0.5 + 0.5
      if (is.list(H)) Hj = H[[j]]
      else Hj = H
      pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
      pred = predict(pi_logit, Hj, s="lambda.min")
      pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
    }
  }
  if(estimate_pi & K==1) {
    Y = (if (is.list(AA)) AA[[1]] else AA)*0.5 + 0.5
    Hj = if (is.list(H)) H[[1]] else H
    pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
    pred = predict(pi_logit, Hj, s="lambda.min")
    pi = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
  }
  
  if (K==1) {
    if (!is.list(AA)) AA = list(AA)
    if (!is.list(RR)) RR = list(RR)
    if (!is.list(pi)) pi = list(pi)
  }
  owl_fit(H, AA, RR, n, K, pi, method, pentype, augment, c, sigma, s, m, solver,
    estimate_pi=estimate_pi)
}

# The same backward recursion also fits training folds and replays validation folds.
owl_fit <- function(H, AA, RR, n, K, pi, method, pentype, augment, c, sigma, s, m, solver,
                    estimate_pi=FALSE, cv=NULL) {
  select=rep(TRUE,n)
  R_future=rep(0,n)
  prob=rep(1,n)
  
  if (augment == T) {
    results = owl_aug(X=H, AA, RR, n, K, pi, method, pentype=pentype, clinear=c, sigma=sigma, s=s, m=m,
      solver=solver, estimate_pi=estimate_pi, cv=cv)
    return(results)
  }
  results=list()
  has_error = 0
  
  for (j in K:1) {
    R = (RR[[j]]+R_future)
    prob = prob*pi[[j]]
    if (method %in% c('hingelinear','hingerbf')) {
      Hj = if (is.list(H)) H[[j]] else H
      if (!is.null(cv) && j==cv$target)
        return(list(H=Hj[select,,drop=FALSE], A=AA[[j]][select], R2=R[select],
          pi=prob[select], models=results))
      if (!is.null(cv$fitted)) {
        results[[j]] = cv$fitted$models[[j]]
        X = Hj[select,,drop=FALSE]
        predict_stage = if (method=='hingelinear') predict.linearcl else predict.rbfcl
        results[[j]]$treatment = if (nrow(X))
          predict_stage(results[[j]], if (nrow(X)==1) as.numeric(X) else X) else numeric(0)
      } else {
        cv.data = if (is.null(cv)) owl_cv_data(H, AA, RR, pi, n, K, j, method, pentype,
          c, sigma, s, m, solver, results, FALSE, estimate_pi) else NULL
        stage_c = if (is.null(cv)) c else cv$tuned[[j]]$c
        stage_sigma = if (is.null(cv)) sigma else cv$tuned[[j]]$sigma
        stage_n = if (is.null(cv)) sum(select) else length(cv$tuned[[j]]$fit)
      }
    }
    
    if (!is.list(H) && is.null(cv$fitted)) {
      if (method=='hingelinear'){
        results[[j]]=owl_single(H[select,,drop=FALSE],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='linear',clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      } else if (method=='hingerbf'){
        results[[j]]=owl_single(H[select,,drop=FALSE],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='rbf',sigma=stage_sigma,clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      } else if (method=='logitlasso'){
        results[[j]]=owl_logit_single(H[select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method, m=m)
      } else if (method=='logit') {
        results[[j]]=owl_logit_single(H[select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method, m=m)
      } else if (method=='ramplinear'){
        results[[j]]=owl_ramp_cv(H[select,], AA[[j]][select], R[select], prob[select], pentype=pentype, bigC=c, bigS=s, K=m)
      } else if (method=='ols' | method=='olslasso') {
        results[[j]]=owl_l2_single(H[select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method, m=m)
      }
    }
    if (is.list(H) && is.null(cv$fitted)) {
      if (method=='hingelinear'){
        results[[j]]=owl_single(H[[j]][select,,drop=FALSE],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='linear',clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      } else if (method=='hingerbf'){
        results[[j]]=owl_single(H[[j]][select,,drop=FALSE],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='rbf',sigma=stage_sigma,clinear=stage_c,m=m,solver=solver,
          cv.data=cv.data,tune=is.null(cv),n_total=stage_n)
      } else if (method=='logitlasso'){
        results[[j]]=owl_logit_single(H[[j]][select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method)
      } else if (method=='logit') {
        results[[j]]=owl_logit_single(H[[j]][select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method)
      } else if (method=='ramplinear'){
        results[[j]]=owl_ramp_cv(H[[j]][select,], AA[[j]][select], R[select], prob[select], pentype=pentype, bigC=c, bigS=s, K=m)
      } else if (method=='ols' | method=='olslasso') {
        results[[j]]=owl_l2_single(H[[j]][select,],AA[[j]][select],R[select],prob[select],pentype=pentype, method=method, m=m)
      }
    }
    if(is.na(max(results[[j]]$fit))) {
      has_error = 1
      break
    }
    select[which(select==TRUE)]=(results[[j]]$treatment==AA[[j]][select])
    
    if ( is.null(cv$fitted) && j>1 && (sum(select)<3 || min(AA[[j-1]][select])==max(AA[[j-1]][select])) ) {
      print(c("Too few observations selected for estimation after stage", j))
      has_error = 1
      break
    }
    R_future=R
  }
  if (method %in% c('hingelinear', 'aughingelinear', 'ramplinear'))  class(results) = 'owl_svmlinear'
  if (method %in% c('hingerbf', 'aughingerbf', 'ramprbf'))  class(results) = 'owl_svmrbf'
  if (method %in% c('logit', 'logitlasso'))  class(results) = 'owl_logit'
  if (method %in% c('ols', 'olslasso'))  class(results) = 'owl_l2'
  if (has_error==1) class(results) = 'has_error'
  class = class(results)
  res = predict(results, H, AA, RR, K, pi)
  return = c(stage=results, list(valuefun=res$valuefun, benefit=res$benefit, pi=pi, type=class))
  class(return) = "owl"
  return
}







