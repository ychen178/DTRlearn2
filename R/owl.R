
#-----------------------------------------------------------------------------------#
# functions for the single-stage (private) and general K-stage outcome-weighted learning
# Yuan Chen, April 2020
#-----------------------------------------------------------------------------------#

## solve single stage outcome-weighted learning
owl_single <-function(H, A, R2, pi, pentype='lasso', kernel='linear', sigma=c(0.03,0.05,0.07), clinear=2.^(-2:2), m=3, e=1e-5)  {
  
  npar=length(clinear)
  n=length(A)
  p=dim(H)[2]
  clinear_input = clinear
  clinear = clinear/n
  
  if (max(R2)!=min(R2)){
    if (pentype=='lasso'){
      cvfit=cv.glmnet(H,R2,nfolds=m)
      co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
    }
    if (pentype=='LSE')  co=coef(lm(R2~H))
    r=R2-cbind(rep(1,n),H)%*%co
  } else r=R2
  r=r/pi
  
  # create CV training and test set
  rand = sample(m,n,replace=T)
  rand_n = 1
  while (rand_n < 10) {
    folds_nobs = rep(0, m)
    for (t in 1:m)  folds_nobs[t] = sum(rand==t)
    folds_nobs = sort(folds_nobs)
    
    if(folds_nobs[1]< 2 | sum(folds_nobs[1:(m-1)]< 3))   rand = sample(1:m, n, replace=T)
    else  break
    rand_n = rand_n + 1
  }
  if (rand_n == 10) {
    print("Too few observations in the training or test set for cross validation. Please increase sample size.")
    return(list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, c=NA, sigma=NA, H=NA, alpha1=NA))
  }
  
  if (kernel=='linear'){
    V=matrix(NA,m,npar)
    for (i in 1:m) {
      this=(rand!=i)
      X=H[this,]
      Y=A[this]
      R=r[this]
      Xt=H[!this,]
      Yt=A[!this]
      Rt=r[!this]
      R_true_t = R2[!this]
      pt=pi[!this]
      for (j in 1:npar) {
        model = wsvm_solve(X, Y, R, C=clinear[j], e=e)
        if (! is.na(max(model$alpha1)) ) {
          YP = predict(model, Xt)
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
      model=wsvm_solve(H,A,r,C=cbest,e=e)
      result = c(model, c=clinear_input[best])
    }
  }
  
  if (kernel=='rbf') {
    nsig=length(sigma)
    V=array(0,c(npar,nsig,m))
    for (i in 1:m){
      this=(rand!=i)
      X=H[this,]
      Y=A[this]
      R=r[this]
      Xt=H[!this,]
      Yt=A[!this]
      Rt=r[!this]
      R_true_t = R2[!this]
      pt=pi[!this]
      for (j in 1:npar){
        for (s in 1:nsig){
          model=wsvm_solve(X,Y,R,'rbf',sigma=sigma[s],C=clinear[j],e=e)
          if (! is.na(max(model$alpha1)) ) {
            YP=predict(model,Xt)
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
      model = wsvm_solve(H,A,r,'rbf', bestSig, C=bestC, e=e)
      result = c(model, c=clinear_input[best[1,1]], sigma=bestSig)
    }
  }
  result
}



## a general K-stage outcome-weighted learning (changed default to augment=FALSE)
owl <- function(H,AA,RR,n,K,pi='estimated', res.lasso=T, loss='hinge', kernel='linear', augment=FALSE, c=2.^(-2:2), sigma=c(0.03,0.05,0.07), s=2.^(-2:2), m=4) {
  
  if(res.lasso==T) pentype = 'lasso'
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
  
  if(pi[[1]][1]=='estimated' & K > 1) {
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
  if(pi[[1]][1]=='estimated' & K==1) {
    Y = AA*0.5 + 0.5
    pi_logit = cv.glmnet(H, Y, family = "binomial", nfolds=3)
    pred = predict(pi_logit, H, s="lambda.min")
    pi = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
  }
  
  if (K==1) {
    if (!is.list(AA)) AA = list(AA)
    if (!is.list(RR)) RR = list(RR)
    if (!is.list(pi)) pi = list(pi)
  }
  select=rep(TRUE,n)
  R_future=rep(0,n)
  prob=rep(1,n)
  
  if (augment == T) {
    results = owl_aug(X=H, AA, RR, n, K, pi, method, pentype=pentype, clinear=c, sigma=sigma, s=s, m=m)
    return(results)
  }
  results=list()
  has_error = 0
  
  for (j in K:1) {
    R = (RR[[j]]+R_future)
    prob = prob*pi[[j]]
    
    if (!is.list(H)) {
      if (method=='hingelinear'){
        results[[j]]=owl_single(H[select,],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='linear',clinear=c,m=m)
      } else if (method=='hingerbf'){
        results[[j]]=owl_single(H[select,],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='rbf',sigma=sigma,clinear=c,m=m)
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
    if (is.list(H)) {
      if (method=='hingelinear'){
        results[[j]]=owl_single(H[[j]][select,],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='linear',clinear=c,m=m)
      } else if (method=='hingerbf'){
        results[[j]]=owl_single(H[[j]][select,],AA[[j]][select],R[select],prob[select],pentype=pentype,kernel='rbf',sigma=sigma,clinear=c,m=m)
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
    
    if ( j>1 && (sum(select)<3 || min(AA[[j-1]][select])==max(AA[[j-1]][select])) ) {
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







