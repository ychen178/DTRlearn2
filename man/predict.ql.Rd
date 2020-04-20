\name{predict.ql}
\alias{predict.ql}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict from a Fitted "ql" Object
}
\description{
This function serves two purposes from a fitted "ql" object. It can recommend treatments for a new independent sample with partial or full subject features observed up to a certain stage. If subject features, treatment assignments and outcomes are fully observed in the new sample, this function can also evaluate the fitted DTR on this new sample, returning the empirical value function and benefit function.
}
\usage{
predict.ql(object, H, AA=NULL, RR=NULL, K, pi=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{fitted "ql" object}

  \item{H}{subject history information before treatment at the \code{K} stages for all subjects in the new sample. It should be constructed the same way as the H in fitting the ql \code{object}. See \code{\link{ql}} for how to construct H. Partial history information is allowed - when first j (j<=K) stages of H is specified, the first j stage treatments will be recommended.}

  \item{AA}{observed treatment assignments at the \code{K} stages for all subjects in the new sample. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. If not specified, treatments will be recommended for the new sample instead of DTR evaluation. The default is \code{AA=NULL}. }

  \item{RR}{observed outcomes at the \code{K} stages for all subjects in the new sample. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. If not specified, treatments will be recommended for the new sample instead of DTR evaluation. The default is \code{RR=NULL}. }

  \item{K}{number of stages of \code{H} observed in the new sample}

  \item{pi}{treatment assignment probabilities of the observed treatments at the K stages for all subjects in the new sample. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. It can be unspecified if one is only interested in treatment recommendations for the new sample. If both \code{AA} and \code{RR} are specified while \code{pi} is not specified, we will estimate the treatment assignment probabilities based on lasso-penalized logistic regressions with predictors being \eqn{H_k} at each stage k. The default is \code{pi=NULL}.}
  
  \item{Qfun}{whether to output the predicted optimal Q-function. The default is Qfun=FALSE. }
  
  \item{fitted}{whether to output the estimated Q-function under the observed treatments. The default is fitted=FALSE. }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
  \item{treatment}{recommennded optimal treatments at the K stages for each subject in the new sample. It is a list of K vectors. If no tailoring variables are learned in the "ql" object, treatments will be assigned randomly with equal probability.}

  \item{valuefun }{overall empirical value function under the fitted DTR evaluated on the new sample. It is returned only when \code{AA} and \code{RR} are fully specified for the K stages.}

  \item{benefit }{overall empirical benefit function under the estimated DTR evaluated on the new sample. It is returned only when \code{AA} and \code{RR} are fully specified for the K stages.}

  \item{pi}{treatment assignment probabilities of the assigned treatments at the K stages for each subject in the new sample. If \code{pi} is not specified but \code{H} and \code{AA} are specified for the K stages, the estimated treatment assignment probabilities from lasso-penalized logistic regressions with \eqn{H_k} being the predictors at each stage k will be returned.}
  
  \item{Q}{the predicted optimal Q-function if Qfun==TRUE.}
  
  \item{fitted}{the estimated Q-function under the observed treatment if fitted==TRUE.}
}
\references{

}
\author{
Yuan Chen, Ying Liu, Donglin Zeng, Yuanjia Wang

Maintainer: Yuan Chen <yc3281@columbia.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ql}}, \code{\link{sim_Kstage}}, \code{\link{owl}}
}
\examples{
# simulate 2-stage training and test sets
n_train = 100
n_test = 500
n_cluster = 10
pinfo = 10
pnoise = 20

train = sim_Kstage(n_train, n_cluster, pinfo, pnoise, K=2)
H1_train = scale(train$X)
H2_train = scale(cbind(H1_train, train$A[[1]], H1_train * train$A[[1]]))
pi_train = list(rep(0.5, n_train), rep(0.5, n_train))

test = sim_Kstage(n_test, n_cluster, pinfo, pnoise, train$centroids, K=2)
H1_test = scale(test$X)
H2_test = scale(cbind(H1_test, test$A[[1]], H1_test * train$A[[1]]))
pi_test = list(rep(0.5, n_test), rep(0.5, n_test))

# estimate DTR with ql on the training sample
ql_train = ql(H=list(H1_train, H2_train), AA=train$A, RR=train$R, K=2, pi=pi_train, m=3)

# evaluate the DTR when full information are observed on the new sample
ql_test = predict(ql_train, H=list(H1_test, H2_test), AA=test$A, RR=test$R, K=2)
ql_test$treatment
ql_test$valuefun
ql_test$pi

# recommned the first-stage treatments only
ql_test2 = predict(ql_train, H=H1_test, K=1)
ql_test2$treatment

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


