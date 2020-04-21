
#-----------------------------------------------------------------------------------#
# single-stage outcome-weighted learning with l2 loss (private function)
# Yuan Chen, June 2018
#-----------------------------------------------------------------------------------#

## solve a single-stage owl with l2 loss
owl_l2_single <-function(H, A, R2, pi, pentype='lasso', method='olslasso', m=4)  {

  return = list()
  n=length(A)
  if (max(R2)!=min(R2)){
    if (pentype=='lasso'){
      cvfit = cv.glmnet(H, R2, nfolds=m)
      co = as.matrix(predict(cvfit, s="lambda.min", type="coeff"))
    } else if (pentype=='LSE'){
      co = coef(lm(R2~H))
    }
    r = R2 - cbind(rep(1,n),H) %*% co
  }
  else r = R2

  r = r/pi
  Y = A*sign(r)

  model = tryCatch(cv.glmnet(H, Y, weights=abs(r), nfolds=m), error=function(e) e)

  if ("error" %in% class(model)) {
    print(model)
    return(list(beta=NA, fit=NA, treatment=NA))
  }
  if (method=='olslasso')
    beta = predict(model, s="lambda.min", type="coeff")
  if (method=='ols')
    beta = predict(model, s=0, type="coeff")

  fit = cbind(1, H) %*% beta
  treatment = 2*(fit > 0) - 1
  results = list(beta=beta, fit=fit, treatment=treatment)
}






