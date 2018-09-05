
#-----------------------------------------------------------------------------------#
# private function for solving the weighted svm
# Yuan Chen, June 2018
#-----------------------------------------------------------------------------------#

## solve the weighted svm via quadratic programming
wsvm_solve <-function(X, A, wR, kernel='linear', sigma=0.05, C=1, e=0.00001) {
  wAR = A * wR
  #print(c("C", C, "sigma", sigma))

  if (kernel=='linear') {
    K = X %*% t(X)
    if (is.vector(X)) K = t(X) %*% X
  }
  else if (kernel=='rbf'){
    rbf = rbfdot(sigma = sigma)
    K = kernelMatrix(rbf, X)
  }

  H = K * (wAR %*% t(wAR))
  eigmat = eigen(H)
  maxev = max(eigmat$values)
  H = eigmat$vectors %*% diag(ifelse(eigmat$values<maxev*1e-5, maxev*1e-5, eigmat$values)) %*% t(eigmat$vectors)

  n = length(A)
  solution = tryCatch(ipop(-abs(wR), H, t(A*wR), 0, numeric(n), rep(C,n), 0, maxiter=100), error=function(e) e)

  if ("error" %in% class(solution)) {
    return(list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, sigma=NA, H=NA, alpha1=NA))
  }
  alpha = primal(solution)
  alpha1 = alpha * wR * A

  if (kernel=='linear'){
    w = t(X) %*% alpha1
    fitted = X %*% w
  } else if (kernel=='rbf'){
    fitted = K %*% alpha1
  }
  rm = sign(wR) * A - fitted
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
    model = list(beta0=bias, beta=w, fit=fit, probability=prob, treatment=2*(fit>0)-1, alpha1=alpha1) #added treatment
    class(model)<-'linearcl'
    } else if (kernel=='rbf') {
    model = list(beta0=bias, fit=fit, probability=prob, treatment=2*(fit>0)-1, sigma=sigma, H=X, alpha1=alpha1) #added treatment, probability
    class(model) = 'rbfcl'
    }
  return (model)
}










