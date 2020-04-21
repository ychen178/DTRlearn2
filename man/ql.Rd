\name{ql}
\alias{ql}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Q-learning for Estimating DTRs
}
\description{
This function implements Q-learning for estimating general K-stage DTRs. Lasso penalty can be applied for variable selection at each stage.
}
\usage{
ql(H, AA, RR, K, pi='estimated', lasso=TRUE, m=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{subject history information before treatment for all subjects at the \code{K} stages. It can be a vector or a matrix when only baseline information is used in estimating the DTR; otherwise, it would be a list of length \code{K}. Please standardize all the variables in \code{H} to have mean 0 and standard deviation 1 before using \code{H} as the input. See details for how to construct H.}

  \item{AA}{observed treatment assignments for all subjects at the \code{K} stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. }

  \item{RR}{observed reward outcomes for all subjects at the \code{K} stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. }

  \item{K}{number of stages}

  \item{pi}{treatment assignment probabilities of the observed treatments for all subjects at the K stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. It can be a user specified input if the treatment assignment probabilities are known.  The default is \code{pi="estimated"}, that is we estimate the treatment assignment probabilities based on lasso-penalized logistic regressions with \eqn{H_k} being the predictors at each stage k.}

  \item{lasso}{specifies whether to add lasso penalty at each stage when fitting the model. The default is \code{lasso=TRUE}. }

  \item{m}{number of folds in the \code{m}-fold cross validation. It is used when \code{res.lasso=T} is specified. The default is \code{m=4}. }

}

\details{
A patient's history information prior to the treatment at stage k can be constructed recursively as \eqn{H_k = (H_{k-1},  A_{k-1}, R_{k-1}, X_k)} with \eqn{H_1=X_1}, where \eqn{X_k} is subject-specific variables collected at stage k just prior to the treatment, \eqn{A_k} is the treatment at stage \eqn{k}, and \eqn{R_k} is the outcome observed post the treatment at stage \eqn{k}. Higher order or interaction terms can also be easily incorporated in \eqn{H_k}, e.g., \eqn{H_k = (H_{k-1}, A_{k-1},R_{k-1}, X_k, H_{k-1}A_{k-1}, R_{k-1}A_{k-1},  X_kA_{k-1})}.
}

\value{
A list of results is returned as an object. It contains the following attributes:
  \item{stage1}{a list of stage 1 results, ...}
  \item{stageK}{a list of stage K results}
  \item{valuefun }{overall empirical value function under the estimated DTR}
  \item{benefit }{overall empirical benefit function under the estimated DTR}
  \item{pi}{treatment assignment probabilities of the assigned treatments for each subject at the K stages. If \code{pi='estimated'} is specified as input, the estimated treatment assignment probabilities from lasso-penalized logistic regressions will be returned. }

In each stage's result, a list is returned which consists of
  \item{co}{the estimated coefficients of \eqn{(1, H, A, H*A)}, the variables in the model at this stage}
  \item{treatment}{the estimated optimal treatment at this stage for each subject in the sample. If no tailoring variables are selected under lasso penalty, treatment will be assigned randomly with equal probability.}
  \item{Q}{the estimated optimal outcome increment from this stage to the end (the estimated optimal Q-function at this stage) for each subject in the sample }
}



\references{
Watkins, C. J. C. H. (1989). Learning from delayed rewards (Doctoral dissertation, University of Cambridge).

Qian, M., & Murphy, S. A. (2011). Performance guarantees for individualized treatment rules. Annals of statistics, 39(2), 1180.
}
\author{
Yuan Chen, Ying Liu, Donglin Zeng, Yuanjia Wang

Maintainer: Yuan Chen <yc3281@columbia.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.ql}}, \code{\link{sim_Kstage}}, \code{\link{owl}}
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

ql_train = ql(H=list(H1_train, H2_train), AA=train$A, RR=train$R, K=2, pi=pi_train, m=3)

ql_test = predict(ql_train, H=list(H1_test, H2_test), AA=test$A, RR=test$R, K=2, pi=pi_test)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
