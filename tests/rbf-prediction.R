# Check batched RBF prediction against independently computed distances.
library(DTRlearn2)
private <- getFromNamespace('predict.rbfcl','DTRlearn2')
set.seed(715)
for (p in c(1L,3L)) {
  training <- matrix(rnorm(12*p),ncol=p)
  queries <- matrix(rnorm(7*p),ncol=p)
  model <- list(beta0=.2, sigma=.3, H=training, alpha1=rnorm(12), svm_orientation=-1)
  class(model) <- 'rbfcl'
  for (rows in list(1L,1:7)) {
    x <- queries[rows,,drop=FALSE]
    kernel <- matrix(NA_real_,nrow(x),nrow(training))
    for (i in seq_len(nrow(x))) for (j in seq_len(nrow(training)))
      kernel[i,j] <- exp(-model$sigma*sum((x[i,]-training[j,])^2))
    expected <- model$beta0+kernel%*%model$alpha1
    treatment <- ifelse(expected>0,1,ifelse(expected<0,-1,1))
    stopifnot(identical(dim(private(model,x)),dim(expected)),
      identical(as.numeric(private(model,x)),as.numeric(treatment)))
    if (nrow(x)==1L)
      stopifnot(identical(as.numeric(private(model,as.numeric(x))),as.numeric(treatment)))
    for (histories in list(x,list(x,x))) {
      object <- list(model,model)
      object$type <- 'owl_svmrbf'
      class(object) <- 'owl'
      for (full_information in c(FALSE,TRUE)) {
        args <- list(object=object,H=histories,K=2)
        if (full_information) {
          args$AA <- rep(list(rep(c(-1,1),length.out=nrow(x))),2)
          args$RR <- rep(list(rep(1,nrow(x))),2)
          args$pi <- rep(list(rep(.5,nrow(x))),2)
        }
        result <- do.call(predict,args)
        for (k in 1:2) {
          stopifnot(identical(dim(result$fit[[k]]),dim(expected)),
            max(abs(result$fit[[k]]-expected))<1e-12,
            identical(as.numeric(result$treatment[[k]]),as.numeric(treatment)))
        }
        if (full_information) {
          matched <- as.numeric(treatment)==args$AA[[1]]
          value <- sum(2*matched/.25)/sum(matched/.25)
          stopifnot(isTRUE(all.equal(result$valuefun,value)))
        }
      }
    }
  }
}
cat('RBF predictions agree with scalar-distance reference, including one-row inputs.\n')
