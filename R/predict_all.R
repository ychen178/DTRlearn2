
#-----------------------------------------------------------------------------------#
# All prediction functions
# Yuan Chen, April 2020
#-----------------------------------------------------------------------------------#

###  prediction functions for intermediate steps (private functions)

predict.linearcl<-function(object,x){
  predict=sign(object$beta0+x%*%object$beta)
  predict
}

predict.rbfcl<-function(object,x){
  rbf=rbfdot(sigma=object$sigma)
  n=dim(object$H)[1]
  if (is.matrix(x)) xm=dim(x)[1]
  else if (is.vector(x) & !is.list(x)) xm=1
  else stop('x must be vector or matrix')
  if (xm==1){ K <- apply(object$H,1,rbf,y=x)
  }else{   K<- matrix(0, xm, n)
  for (i in 1:xm) {K[i,]=apply(object$H,1,rbf,y=x[i,]) }}
  predict=sign(object$beta0+K%*%object$alpha1)
  predict
}


###  predict with owl() objects
predict.owl <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL, ...) {
  
  # estimate pi if AA have inputs
  if(!is.null(AA) & is.null(pi)) {
    pi = list()
    for (j in 1:K) {
      Y = AA[[j]]*0.5 + 0.5
      if (is.list(H)) Hj = H[[j]]
      pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
      pred = predict(pi_logit, Hj, s="lambda.min")
      pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
    }
  }
  
  res = list()
  if(object$type == "owl_svmlinear") res = predict.owl_svmlinear(object, H, AA, RR, K, pi)
  else if(object$type == "owl_svmrbf") res = predict.owl_svmrbf(object, H, AA, RR, K, pi)
  else if(object$type == "owl_logit") res = predict.owl_logit(object, H, AA, RR, K, pi)
  else if(object$type == "owl_l2") res = predict.owl_l2(object, H, AA, RR, K, pi)
  else if(object$type == "has_error") res = predict.has_error(object)
  res
}


###  predict with  hingelinear, aughingelinear object
predict.owl_svmlinear <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL) {
  
  if (is.vector(H) & !is.list(H))  H = matrix(H, ncol=1)
  
  fit = list()
  predprob = list()
  treatment = list()
  
  if (is.null(AA) | is.null(RR)) { #only H and K available
    for (i in 1:K) {
      if (is.matrix(H))
        fit[[i]] = object[[i]]$beta0 + H %*% object[[i]]$beta
      else if (is.list(H))
        fit[[i]] = object[[i]]$beta0 + H[[i]] %*% object[[i]]$beta
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      treatment[[i]] = sign(fit[[i]])
      predprob[[i]] = exp(fit[[i]]) / (1+ exp(fit[[i]]))
    }
    return = list(fit = fit, probability=predprob, treatment = treatment)
  }
  else {  # full information sample
    if(is.null(pi)) {  # estimate pi
      pi = list()
      for (j in 1:K) {
        Y = AA[[j]]*0.5 + 0.5
        if (is.list(H)) Hj = H[[j]]
        pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
        pred = predict(pi_logit, Hj, s="lambda.min")
        pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
      }
    }
    if (K==1) {
      if (!is.list(AA)) AA = list(AA)
      if (!is.list(RR)) RR = list(RR)
      if (!is.list(pi)) pi = list(pi)
    }
    n = length(AA[[1]])
    select = rep(1,n)
    prob = rep(1,n)
    sumR = rep(0,n)
    
    for (i in 1:K) {
      if (is.matrix(H))
        fit[[i]] = object[[i]]$beta0 + H %*% object[[i]]$beta
      else if (is.list(H))
        fit[[i]] = object[[i]]$beta0 + H[[i]] %*% object[[i]]$beta
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      
      treatment[[i]] = sign(fit[[i]])
      predprob[[i]] = exp(fit[[i]]) / (1+ exp(fit[[i]]))
      select = select * (treatment[[i]] == AA[[i]])
      prob = prob * pi[[i]]
      sumR = sumR + RR[[i]]
    }
    valuefun = sum(sumR * select / prob) / n
    benefitfun = valuefun - sum(sumR*(1-select)/prob)/n
    return = list(fit = fit, treatment = treatment, probability=predprob, valuefun = valuefun, benefit = benefitfun, pi=pi)
  }
  return
}

###  predict with  hingerbf, aughingerbf   object
predict.owl_svmrbf <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL) {
  
  if (is.vector(H) & !is.list(H))  H = matrix(H, ncol=1)
  fit = list()
  predprob = list()
  treatment = list()
  
  if (is.null(AA) | is.null(RR)) {
    for(i in 1:K) {
      rbf = rbfdot(sigma = object[[i]]$sigma)
      nsv = dim(object[[i]]$H)[1]
      if (is.matrix(H))   n = dim(H)[1]
      else if (is.list(H))   n = dim(H[[i]])[1]
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      
      kern = matrix(0, n, nsv)
      for (j in 1 : n) {
        if (is.matrix(H))
          kern[j,] = apply(object[[i]]$H, 1, rbf, y=H[j,])
        else if (is.list(H))
          kern[j,] = apply(object[[i]]$H, 1, rbf, y=H[[i]][j,])
      }
      fit[[i]] = object[[i]]$beta0 + kern %*% object[[i]]$alpha1
      predprob[[i]] = exp(fit[[i]]) / (1+ exp(fit[[i]]))
      treatment[[i]] = sign(fit[[i]])
    }
    return = list(fit = fit, probability=predprob, treatment = treatment)
  }
  else {  # full information sample
    if(is.null(pi)) {  # estimate pi
      pi = list()
      for (j in 1:K) {
        Y = AA[[j]]*0.5 + 0.5
        if (is.list(H)) Hj = H[[j]]
        pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
        pred = predict(pi_logit, Hj, s="lambda.min")
        pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
      }
    }
    if (K==1) {
      if (!is.list(AA)) AA = list(AA)
      if (!is.list(RR)) RR = list(RR)
      if (!is.list(pi)) pi = list(pi)
    }
    n = length(AA[[1]])
    select = rep(1,n)
    prob = rep(1,n)
    sumR = rep(0, n)
    
    for(i in 1:K) {
      rbf = rbfdot(sigma = object[[i]]$sigma)
      if(is.matrix(object[[i]]$H))  nsv = dim(object[[i]]$H)[1]
      if(is.vector(object[[i]]$H))  nsv = length(object[[i]]$H)
      kern = matrix(0, n, nsv)
      for (j in 1 : n) {
        if (is.matrix(H))
          kern[j,] = apply(object[[i]]$H, 1, rbf, y=H[j,])
        else if (is.list(H))
          kern[j,] = apply(object[[i]]$H, 1, rbf, y=H[[i]][j,])
        else stop(gettextf("H must be a vector or matrix, or a list of vectors or matrice"))
      }
      fit[[i]] = object[[i]]$beta0 + kern %*% object[[i]]$alpha1
      predprob[[i]] = exp(fit[[i]]) / (1+exp(fit[[i]]))
      treatment[[i]] = sign(fit[[i]])
      select = select * (treatment[[i]] == AA[[i]])
      prob = prob * pi[[i]]
      sumR = sumR + RR[[i]]
    }
    valuefun = sum(sumR * select / prob) / n
    benefitfun = valuefun - sum(sumR*(1-select)/prob)/n
    return = list(fit = fit, treatment = treatment, probability=predprob, valuefun = valuefun, benefit = benefitfun, pi=pi)
  }
  return
}



### predict with  logit / logitlasso  object
predict.owl_logit <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL) {
  
  if (is.vector(H) & !is.list(H))   H = matrix(H, ncol=1)
  fit = list()
  predprob = list()
  treatment = list()
  
  if (is.null(AA) | is.null(RR)) {
    for (i in 1:K) {
      if (is.matrix(H)) {
        n = dim(H)[1]
        xbeta = cbind(rep(1,n), H) %*% object[[i]]$beta
      } else if (is.list(H)) {
        if(is.matrix(H[[i]]))  n = dim(H[[i]])[1]
        else if (is.vector(H[[i]]))  n = length(H[[i]])
        xbeta = cbind(rep(1,n), H[[i]]) %*% object[[i]]$beta
      }
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      fit[[i]] = xbeta
      predprob[[i]] = exp(xbeta)/(1 + exp(xbeta))
      treatment[[i]] = 2 * (predprob[[i]] > 0.5) - 1
    }
    return = list(fit=fit, probability=predprob, treatment = treatment)
  }
  else {  # full information sample
    if(is.null(pi)) {  # estimate pi
      pi = list()
      for (j in 1:K) {
        Y = AA[[j]]*0.5 + 0.5
        if (is.list(H)) Hj = H[[j]]
        pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
        pred = predict(pi_logit, Hj, s="lambda.min")
        pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
      }
    }
    if (K==1) {
      if (!is.list(AA)) AA = list(AA)
      if (!is.list(RR)) RR = list(RR)
      if (!is.list(pi)) pi = list(pi)
    }
    n = length(AA[[1]])
    select = rep(1,n)
    prob = rep(1,n)
    sumR = rep(0, n)
    
    for (i in 1:K) {
      if (is.matrix(H))
        xbeta = cbind(rep(1,n), H) %*% object[[i]]$beta
      else if (is.list(H))
        xbeta = cbind(rep(1,n), H[[i]]) %*% object[[i]]$beta
      else stop(gettextf("H must be a vector or matrix, or a list of vectors or matrice"))
      
      fit[[i]] = xbeta
      predprob[[i]] = exp(xbeta)/(1 + exp(xbeta))
      treatment[[i]] = 2 * (predprob[[i]] > 0.5) - 1
      select = select * (treatment[[i]] == AA[[i]])
      prob = prob * pi[[i]]
      sumR = sumR + RR[[i]]
    }
    valuefun = sum(sumR * select / prob, na.rm = T) / sum(!is.na(sumR * select)) 
    benefitfun = valuefun - sum(sumR*(1-select)/prob, na.rm = T) / sum(!is.na(sumR * select)) 
    return = list(fit = xbeta, probability =predprob, treatment = treatment, valuefun = valuefun, benefit = benefitfun, pi=pi)
  }
  return
}


### predict based on   l2 / l2.lasso  owl() object
predict.owl_l2 <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL) {
  
  if (is.vector(H) & !is.list(H))   H = matrix(H, ncol=1)
  fit = list()
  predprob = list()
  treatment = list()
  
  if (is.null(AA) | is.null(RR)) {
    for (i in 1:K) {
      if (is.matrix(H)) {
        n = dim(H)[1]
        fit[[i]] = cbind(rep(1,n), H) %*% object[[i]]$beta
      }
      else if (is.list(H)) {
        if(is.matrix(H[[i]]))  n = dim(H[[i]])[1]
        else if (is.vector(H[[i]]))  n = length(H[[i]])
        fit[[i]] = cbind(rep(1,n), H[[i]]) %*% object[[i]]$beta
      }
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      predprob[[i]] = exp(fit[[i]])/(1 + exp(fit[[i]]))
      treatment[[i]] = 2*(fit[[i]] > 0) - 1
    }
    return = list(fit = fit, probability=predprob, treatment = treatment)
  }
  else {   # full information sample
    if(is.null(pi)) {  # estimate pi
      pi = list()
      for (j in 1:K) {
        Y = AA[[j]]*0.5 + 0.5
        if (is.list(H)) Hj = H[[j]]
        pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
        pred = predict(pi_logit, Hj, s="lambda.min")
        pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
      }
    }
    if (K==1) {
      if (!is.list(AA)) AA = list(AA)
      if (!is.list(RR)) RR = list(RR)
      if (!is.list(pi)) pi = list(pi)
    }
    n = length(AA[[1]])
    select = rep(1,n)
    prob = rep(1,n)
    sumR = rep(0,n)
    
    for (i in 1:K) {
      if (is.matrix(H))
        fit[[i]] = cbind(rep(1,n), H) %*% object[[i]]$beta
      else if (is.list(H))
        fit[[i]] = cbind(rep(1,n), H[[i]]) %*% object[[i]]$beta
      else stop(gettextf("H must be a vector or matrix, or a list of vectors or matrice"))
      
      predprob[[i]] = exp(fit[[i]])/(1 + exp(fit[[i]]))
      treatment[[i]] = 2*(fit[[i]] > 0) - 1
      select = select * (treatment[[i]] == AA[[i]])
      prob = prob * pi[[i]]
      sumR = sumR + RR[[i]]
    }
    valuefun = sum(sumR * select / prob) / n
    benefitfun = valuefun - sum(sumR*(1-select)/prob)/n
    return = list(fit = fit, probability=predprob, treatment = treatment, valuefun = valuefun, benefit = benefitfun, pi=pi)
  }
  return
}

## handle the NA case (when there is error in the fitted object)
predict.has_error = function(object) {
  list(fit = NA, treatment = NA, valuefun = NA, benefit = NA)
}


### predict with ql() object
predict.ql <- function(object, H, AA=NULL, RR=NULL, K, pi=NULL, Qopt=FALSE, Qfit=FALSE, ...) {  # added whether to predict the fitted Q function under the given treatment
  
  if (is.vector(H) & !is.list(H))   H = matrix(H, ncol=1)
  Q = list()      # predicted optimal Q function
  fit = list()
  treatment = list()
  fitted = list()  # fitted Q function
  
  # estimate pi if AA have inputs
  if(!is.null(AA) & is.null(pi)) {
    pi = list()
    for (j in 1:K) {
      Y = AA[[j]]*0.5 + 0.5
      if (is.list(H)) Hj = H[[j]]
      pi_logit = cv.glmnet(Hj, Y, family = "binomial", nfolds=3)
      pred = predict(pi_logit, Hj, s="lambda.min")
      pi[[j]] = exp(pred)/(1 + exp(pred)) * (Y==1) + 1/(1 + exp(pred)) * (Y==0)
    }
  }

  if (K==1) {
    if (!is.null(AA) & !is.list(AA)) AA = list(AA)
    if (!is.null(RR) & !is.list(RR)) RR = list(RR)
    if (!is.null(pi) & !is.list(pi)) pi = list(pi)
  }
  
  if (is.null(RR)) { #no outcome
    for (i in 1:K) {
      if (is.matrix(H)) {
        p = dim(H)[2]
        n = dim(H)[1]
        fit[[i]] = (object[[i]]$co[p+2] + H %*% object[[i]]$co[(p+3):(2*p+2)])
        treatment[[i]] = 2*(fit[[i]]>0) - 1
        if(min(fit[[i]])==max(fit[[i]])) treatment[[i]] = rbinom(n, 1, 0.5) #no tailoring vars, randomize treatments
        if(Qopt==TRUE)  Q[[i]] = cbind(rep(1,n), H, c(as.vector(treatment[[i]])), diag(c(as.vector(treatment[[i]]))) %*% H) %*% object[[i]]$co
        if(Qfit==TRUE & !is.null(AA)) fitted[[i]] = cbind(rep(1,n), H, c(as.vector(AA[[i]])), diag(c(as.vector(AA[[i]]))) %*% H) %*% object[[i]]$co
      }
      else if (is.list(H)) {
        p = dim(H[[i]])[2]
        n = dim(H[[i]])[1]
        fit[[i]] = (object[[i]]$co[p+2] + H[[i]] %*% object[[i]]$co[(p+3):(2*p+2)])
        treatment[[i]] = 2*(fit[[i]]>0) - 1
        if(min(fit[[i]])==max(fit[[i]])) treatment[[i]] = rbinom(n, 1, 0.5)
        if(Qopt==TRUE)  Q[[i]] = cbind(rep(1,n), H[[i]], c(as.vector(treatment[[i]])), diag(c(as.vector(treatment[[i]]))) %*% H[[i]]) %*% object[[i]]$co
        if(Qfit==TRUE & !is.null(AA)) fitted[[i]] = cbind(rep(1,n), H[[i]], c(as.vector(AA[[i]])), diag(c(as.vector(AA[[i]]))) %*% H[[i]]) %*% object[[i]]$co 
      }
      else stop(gettextf("H must be a vector or matrix, or a list of vectors or matrices"))
      if(min(fit[[i]])==max(fit[[i]])) treatment[[i]] = rbinom(n, 1, 0.5) #no tailoring vars, randomize treatments
    }
    if(Qopt==TRUE | Qfit==TRUE)  return = list(treatment = treatment, Q=Q, fitted=fitted)
    else return = list(treatment = treatment)
  }
  else {  # with outcome
    n = length(AA[[1]])
    select = rep(1,n)
    prob = rep(1,n)
    sumR = rep(0, n)
    
    for (i in 1:K) {
      if (is.matrix(H)) {
        p = dim(H)[2]
        fit[[i]] = (object[[i]]$co[p+2] + H %*% object[[i]]$co[(p+3):(2*p+2)])
        treatment[[i]] = 2*(fit[[i]]>0) - 1
        if(Qopt==TRUE)  Q[[i]] = cbind(rep(1,n), H, c(as.vector(treatment[[i]])), diag(c(as.vector(treatment[[i]]))) %*% H) %*% object[[i]]$co
        if(Qfit==TRUE & !is.null(AA)) fitted[[i]] = cbind(rep(1,n), H, c(as.vector(AA[[i]])), diag(c(as.vector(AA[[i]]))) %*% H) %*% object[[i]]$co
      } else if (is.list(H)) {
        p = dim(H[[i]])[2]
        fit[[i]] = (object[[i]]$co[p+2] + H[[i]] %*% object[[i]]$co[(p+3):(2*p+2)])
        treatment[[i]] = 2*(fit[[i]]>0) - 1
        if(Qopt==TRUE)  Q[[i]] = cbind(rep(1,n), H[[i]], c(as.vector(treatment[[i]])), diag(c(as.vector(treatment[[i]]))) %*% H[[i]]) %*% object[[i]]$co 
        if(Qfit==TRUE & !is.null(AA)) fitted[[i]] = cbind(rep(1,n), H[[i]], c(as.vector(AA[[i]])), diag(c(as.vector(AA[[i]]))) %*% H[[i]]) %*% object[[i]]$co 
      }
      else stop(gettextf("H must be a vector or matrix, or a list of K vectors or matrice"))
      
      if(min(fit[[i]])==max(fit[[i]])) treatment[[i]] = rbinom(n, 1, 0.5) #no tailoring vars, randomize treatments
      select = select * (treatment[[i]] == AA[[i]])
      prob = prob * pi[[i]]
      sumR = sumR + RR[[i]]
    }
    valuefun = mean(sumR * select / prob)
    benefitfun = valuefun - mean(sumR*(1-select)/prob)
    if(Qopt==TRUE | Qfit==TRUE) return = list(treatment = treatment, valuefun = valuefun, benefit = benefitfun, pi=pi, Q=Q, fitted=fitted)
    else  return = list(treatment = treatment, valuefun = valuefun, benefit = benefitfun, pi=pi)
  }
  return
}




## handle the NA case (when there is error in the fitted object)
predict.has_error = function(object) {
  list(fit = NA, treatment = NA, valuefun = NA, benefit = NA)
}




