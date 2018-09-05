
#-----------------------------------------------------------------------------------#
# Simulate a K-stage SMART data with n_cluster clusters characterized by the differential means in the pinfo variables
# Yuan Chen, June 2018
#-----------------------------------------------------------------------------------#


sim_Kstage = function(n, n_cluster, pinfo, pnoise, centroids=NULL, K) {

  if(is.null(centroids))
    centroids = matrix(rnorm(n_cluster * pinfo, 0, sqrt(5)), nrow = n_cluster)
  y=list()
  R=list()
  A=list()
  X=matrix(0, n, pinfo+pnoise)
  sigma1 = ifelse(diag(pinfo), 1, 0.2) #covariance matrix for Xpinfo

  for (k in 1:K){
    y[[k]] = rep(0,n)
    R[[k]] = rep(0,n)
    A[[k]] = 2*rbinom(n, 1, 0.5) - 1
  }

  # assign optimal treatment y and X, based on cluster
  end = 1
  for (l in 1:n_cluster){
    start = end
    end = end + n/n_cluster
    for (k in 1:K) {
      y[[k]][start:(end-1)] = 2 * (floor(l/(2*k-1)) %% 2) - 1
    }
    X[start:(end-1),1:pinfo] =
      mvrnorm(n/n_cluster, centroids[l,], Sigma = sigma1)
  }
  X[, pinfo + (1:pnoise)] = mvrnorm(n, rep(0,pnoise), diag(pnoise))

  for(k in 1:K) {
    R[[K]] = R[[K]] + A[[k]]*y[[k]]
  }
  R[[K]] = R[[K]] + rnorm(n)

  list(X=X, A=A, R=R, optA=y, centroids=centroids)
}




