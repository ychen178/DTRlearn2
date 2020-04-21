
#-----------------------------------------------------------------------------------#
# single-stage outcome-weighted learning with binomial deviance loss (private function)
# Yuan Chen, June 2018
#-----------------------------------------------------------------------------------#

## solve single-stage owl with logit loss
owl_logit_single <- function(H,A,R2,pi,pentype='lasso',method='logitlasso', m=4) {

  n=length(A)
  if (max(R2)!=min(R2)){
    if (pentype=='lasso'){
      cvfit=cv.glmnet(H,R2,nfolds=m)
      co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
    }else if (pentype=='LSE'){
      co=coef(lm(R2~H))
    }
    r=R2-cbind(rep(1,n),H)%*%co
  }
  else r=R2

  r=r/pi
  Y=A*sign(r)*0.5+0.5

  cvlogit = tryCatch(cv.glmnet(H, Y, family = "binomial", weights=abs(r), nfolds=m), error=function(e) e)
  if("error" %in% class(cvlogit))  {
    has_error = 1
    print(cvlogit)
    beta = NA; xbeta = NA; prob = NA; tx = NA
  } else {
    if (method=='logitlasso')
      beta = predict(cvlogit,s="lambda.min",type="coeff")
    if (method=='logit')
      beta = predict(cvlogit, s=0, type="coeff")

    xbeta = cbind(rep(1,n), H) %*% beta
    prob = exp(xbeta)/(1 + exp(xbeta))
    tx = 1*(prob>=0.5) - 1*(prob<0.5)
  }
  results = list(beta=beta, fit=xbeta, probability=prob, treatment=tx)
}



