\name{owl}
\alias{owl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated Outcome-weighted Learning for Estimating DTRs
}
\description{
This function implements a variety of outcome-weighted learning methods for estimating general K-stage DTRs. Different loss functions - SVM hinge loss, SVM ramp loss, binomial deviance loss, and L2 loss - can be adopted to solve the weighted classification problem at each stage. Augmentation in the outcomes is allowed to improve efficiency especially when there are multiple stages with a small sample size. Cross validation is conducted to choose the best tuning parameters if any.
}
\usage{
owl(H, AA, RR, n, K, pi='estimated', res.lasso=TRUE, loss='hinge', kernel='linear',
augment=FALSE, c=2^(-2:2), sigma=c(0.03,0.05,0.07), s=2.^(-2:2), m=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{subject history information before treatment for the \code{K} stages. It can be a vector or a matrix when only baseline information is used in estimating the DTR; otherwise, it would be a list of length \code{K}. Please standardize all the variables in \code{H} to have mean 0 and standard deviation 1 before using \code{H} as the input. See details for how to construct H.}

  \item{AA}{observed treatment assignments for all subjects at the \code{K} stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. }

  \item{RR}{observed reward outcomes for all subjects at the \code{K} stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. }

  \item{n}{sample size, number of subjects in the dataset}

  \item{K}{number of stages}

  \item{pi}{treatment assignment probabilities of the observed treatments for all subjects at the K stages. It is a vector if \code{K=1}, or a list of \code{K} vectors corresponding to the \code{K} stages. It can be a user specified input if the treatment assignment probabilities are known.  The default is \code{pi="estimated"}, that is we estimate the treatment assignment probabilities based on lasso-penalized logistic regressions with \eqn{H_k} being the predictors at each stage k.}

  \item{res.lasso}{whether or not to use lasso penalty in the regression to take residuals for constructing the weights. The default is \code{res.lasso=TRUE}.}

  \item{loss}{loss function for sovling the weighted classification problem at each stage. The options are \code{"hinge", "ramp", "logit", "logit.lasso", "l2", "l2.lasso"}. \code{"hinge"} and \code{"ramp"} are for the SVM hinge loss and SVM ramp loss. \code{"logit"} and \code{"logit.lasso"} are for the binomial deviance loss used in the logistic regression, where lasso penalty is applied under \code{"logit.lasso"}. \code{"l2"} and \code{"l2.lasso"} are for the L2 or square loss, where lasso penalty is applied under \code{"l2.lasso"}. The default is \code{loss="hinge"}. }

  \item{kernel}{kernel function to use under SVM hinge loss or SVM ramp loss. \code{"linear"} and \code{"rbf"} kernel are implemented under SVM hinge loss; \code{"linear"} kernel is implemented under SVM ramp loss. The default is \code{kernel="linear"}.}

  \item{augment}{whether or not to use augmented outcomes at each stage. Augmentation is recommended when there are multiple stages and the sample size is small. The default is \code{augment=FALSE}.}

  \item{c}{a vector specifies the values of the regularization parameter C for tuning under SVM hinge loss or SVM ramp loss. The default is \code{c=2^(-2:2)}. In practice, a wider range of \code{c} can be specified based on the data. }

  \item{sigma}{a vector specifies the values of the positive parameter \eqn{\sigma}{sigma} in the RBF kernel for tuning under SVM hinge loss, i.e., when \code{loss="hinge"} and \code{kernel="rbf"}. The default is  \code{sigma=c(0.03,0.05,0.07)}. In practice, a wider range of \code{sigma} can be specified based on the data. }

  \item{s}{a vector specifies the values of the slope parameter in the SVM ramp loss for tuning, i.e., when \code{loss="ramp"} and \code{kernel="linear"}.
  The default is \code{c=2^(-2:2)}. In practice, a wider range of \code{s} can be specified based on the data. }

  \item{m}{number of folds in the m-fold cross validation for choosing the tuning parameters \code{c}, \code{sigma} or \code{s}. It is also used for choosing the tuning parameter of the lasso penalty when \code{res.lasso=T}, \code{loss="logit.lasso"} or \code{loss="l2.lasso"} is specified. The default is \code{m=4}.}

}


\details{
A patient's history information prior to the treatment at stage k can be constructed recursively as \eqn{H_k = (H_{k-1},  A_{k-1}, R_{k-1}, X_k)} with \eqn{H_1=X_1}, where \eqn{X_k} is subject-specific variables collected at stage k just prior to the treatment, \eqn{A_k} is the treatment at stage \eqn{k}, and \eqn{R_k} is the outcome observed post the treatment at stage \eqn{k}. Higher order or interaction terms can also be easily incorporated in \eqn{H_k}, e.g., \eqn{H_k = (H_{k-1}, A_{k-1}, R_{k-1}, X_k, H_{k-1}A_{k-1}, R_{k-1}A_{k-1},  X_kA_{k-1})}.
}

\value{
A list of results is returned as an object. It contains the following attributes:
  \item{stage1 }{a list of stage 1 results, ...}
  \item{stageK }{a list of stage K results}
  \item{valuefun }{overall empirical value function under the estimated DTR}
  \item{benefit }{overall empirical benefit function under the estimated DTR}
  \item{pi}{treatment assignment probabilities of the observed treatments for each subject at the K stages. It is a list of K vectors. If \code{pi='estimated'} is specified as input, the estimated treatment assignment probabilities from lasso-penalized logistic regressions will be returned. }
  \item{type}{object type corresponding to the specified \code{loss} and \code{kernel} }


In each stage's result, a list is returned which consists of
  \item{beta0}{estimated coefficient of the intercept in the decision function}
  \item{beta}{estimated coefficients of \eqn{H_k} in the decision function. It's not returned with RBF kernel under SVM hinge loss.}
  \item{fit}{fitted decision function for each subject}
  \item{probability}{estimated probability that treatment 1 (vs. -1) is the optimal treatment for each subject in the sample. It's calculated by exp(fit)/(1 + exp(fit)).  }
  \item{treatment}{the estimated optimal treatment for each subject}
  \item{c}{the best regularization parameter C in SVM hinge loss or SVM ramp loss, chosen from the values specified in \code{c} via cross validation }

  \item{sigma}{the best parameter \eqn{\sigma} in the RBF kernel, chosen from the values specified in \code{sigma} via cross validation }

  \item{s}{the best slope parameter s in the ramp loss, chosen from the values specified in \code{s} via cross validation. }
  \item{iter}{number of iterations under SVM ramp loss}

  \item{alpha1}{the solution to the Lagrangian dual problem under SVM hinge loss or SVM ramp loss. It is used for constructing the decision function on the new sample. }

  \item{H}{the input H, returned only under SVM hinge loss with RBF kernel. It is used for constructing the RBF kernel on the new sample. }
  }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{
Liu, Y., Wang, Y., Kosorok, M., Zhao, Y., & Zeng, D. (2014). Robust hybrid learning for estimating personalized dynamic treatment regimens. arXiv preprint. arXiv, 1611.

Liu, Y., Wang, Y., Kosorok, M., Zhao, Y., & Zeng, D. (2018). Augmented Outcome-weighted Learning for Estimating Optimal Dynamic Treatment Regimens. Statistics in Medicine. In press.

Zhao, Y., Zeng, D., Rush, A. J., & Kosorok, M. R. (2012). Estimating individualized treatment rules using outcome weighted learning. Journal of the American Statistical Association, 107(499), 1106-1118.

Zhao, Y. Q., Zeng, D., Laber, E. B., & Kosorok, M. R. (2015). New statistical learning methods for estimating optimal dynamic treatment regimes. Journal of the American Statistical Association, 110(510), 583-598.

}
\author{
Yuan Chen, Ying Liu, Donglin Zeng, Yuanjia Wang

Maintainer: Yuan Chen <yc3281@columbia.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.owl}}, \code{\link{sim_Kstage}}, \code{\link{ql}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

# estimate DTR with owl on the training sample
owl_train = owl(H=list(H1_train, H2_train), AA=train$A, RR=train$R, n=n_train, K=2, pi=pi_train,
    loss='hinge', augment=TRUE, m=3)
owl_train$stage1$beta
owl_train$stage1$treatment
owl_train$valuefun

# apply the estimated DTR to the test sample
owl_test = predict(owl_train, H=list(H1_test, H2_test), AA=test$A, RR=test$R, K=2, pi=pi_test)
owl_test$treatment
owl_test$valuefun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

