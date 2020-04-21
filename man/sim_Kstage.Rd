\name{sim_Kstage}
\alias{sim_Kstage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a K-stage SMART data
}
\description{
This function simulates a K-stage SMART data with \code{(pinfo + pnoise)} baseline variables from a multivariate Gaussian distribution. The \code{pinfo} variables have variance 1 and pairwise correlation 0.2; the \code{pnoise} variables have mean 0 and are uncorrelated with each other and with the \code{pinfo} variables.

Subjects are from \code{n_cluster} latent groups with equal sizes, and these \code{n_cluster} groups are characterized by their differentiable means in the \code{pinfo} feature variables. Each latent group has its own optimal treatment sequence, where the optimal treatment for subjects in group g at stage k is generated as \eqn{ A^* = 2( [ g/(2k -1) ]} mod 2) - 1. The assigned treatment group (1 or -1) for each subject at each stage is randomly generated with equal probability. The primary outcome is observed only at the end of the trial, which is generated as
\eqn{R = \sum_{k=1}^{K} A_k A_k^* + N(0,1)}.
}

\usage{
sim_Kstage (n, n_cluster, pinfo, pnoise, centroids=NULL, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size, should be a multiple of \code{n_cluster}.}
  \item{n_cluster}{number of latent groups}
  \item{pinfo}{number of informative baseline variables}
  \item{pnoise}{number of non-informative baseline variables}
  \item{centroids}{centroids of the \code{pinfo} variables for the \code{n_cluster} groups. It is a matrix of dimension \code{n_cluster} by \code{pinfo}. It's used as the means of the multivariate Gaussians to generate the \code{pinfo} variables for the \code{n_cluster} groups. For a training set, do not assign centroids, the centroids are generated randomly from N(0,5) by the function. For a test set, one should assign the same set of centroids as the training set.}
  \item{K}{number of stages.}
}

\value{
  \item{X }{baseline variables. It is a matrix of dimension \code{n} by \code{(pinfo + pnoise)}.}
  \item{A }{treatment assigments for the K-stages. It is a list of K vectors.}
  \item{R }{outcomes of the K-stages. It is a list of K vectors. In this simulation setting, no intermediate outcomes are observed, so the first K-1 vectors are vectors of 0.}
  \item{optA }{optimal treatments for the K-stages. It is a list of K vectors.}
  \item{centroids }{centroids of the \code{pinfo} variables for the \code{n_cluster} groups. It is a matrix of dimension \code{n_cluster} by \code{pinfo}. }

}

\author{
Yuan Chen, Ying Liu, Donglin Zeng, Yuanjia Wang

Maintainer: Yuan Chen <yc3281@columbia.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{owl}}, \code{\link{ql}}
}
\examples{
n_train = 100
n_test = 500
n_cluster = 10
pinfo = 10
pnoise = 20

# simulate a 2-stage training set
train = sim_Kstage(n_train, n_cluster, pinfo, pnoise, K=2)

# simulate an independent 2-stage test set with the same centroids of the training set
test = sim_Kstage(n_test, n_cluster, pinfo, pnoise, train$centroids, K=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




