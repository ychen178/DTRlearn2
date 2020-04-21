#-----------------------------------------------------------------------------------#
# private function for solving the weighted svm with quadratic programming
# Yuan Chen, April 2020
#-----------------------------------------------------------------------------------#

wsvm_solve <-function(X, A, wR, kernel='linear', sigma=0.05, C=1, e=1e-7) {
  
  if (kernel=='linear') {
    K = X %*% t(X)
    if (is.vector(X)) K = t(X) %*% X
  }
  else if (kernel=='rbf'){
    rbf = rbfdot(sigma = sigma)
    K = kernelMatrix(rbf, X)
  }
  
  y = A * sign(wR)
  H = y %*% t(y) * K
  H = H + 1e-8 * diag(NCOL(K)) %*% (tcrossprod(wR))
  
  
  n = length(A)
  solution <- tryCatch(ipop(c = rep(-1, n), H = H, A = t(y), b = 0, l = numeric(n), u = C*abs(wR), r = 0), error=function(er) er)
  if ("error" %in% class(solution)) {
    return(list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, sigma=NA, H=NA, alpha1=NA))
  }
  alpha = primal(solution)
  alpha1 = alpha * y 
  
  if (kernel=='linear'){
    w = t(X) %*% alpha1
    fitted = X %*% w
  } else if (kernel=='rbf'){
    fitted = K %*% alpha1
  }
  rm = y - fitted
  Imid = (alpha < C-e) & (alpha > e)
  rmid = rm[Imid==1]
  if (sum(Imid)>0){
    bias = mean(rmid)
  } else {
    Iup = ((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
    Ilow = ((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
    rup = rm[Iup]
    rlow = rm[Ilow]
    bias = (min(rup)+max(rlow))/2
  }
  fit = bias + fitted
  prob = exp(fit) / (1+ exp(fit))
  
  
  if (kernel=='linear') {
    model = list(beta0=bias, beta=w, fit=fit, probability=prob, treatment=2*(fit>0)-1, alpha1=alpha1) #, solution=solution) 
    class(model)<-'linearcl'
  } else if (kernel=='rbf') {
    model = list(beta0=bias, fit=fit, probability=prob, treatment=2*(fit>0)-1, sigma=sigma, H=X, alpha1=alpha1)
    class(model) = 'rbfcl'
  }
  return (model)
}




