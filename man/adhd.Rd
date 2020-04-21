\name{adhd}
\alias{adhd}
\docType{data}
\title{
A 2-stage SMART data of children with ADHD
}
\description{
We provide a two-stage sequential multiple assignment randomized trial (SMART) data of 150 children with ADHD mimicking a real world study. At the first stage, children were randomized to treatment of low-intensity behavioral modification (BMOD) or low-intensity methamphetamine (MED) with equal probability. At second stage, children were randomized to treatment of low-intensity BMOD + low-intensity MED, or high-intensity BMOD with equal probability. The primary outcome of study was children's school performance score ranging from 1 to 5 assessed at the end of the study for all participants.
}
\usage{data("adhd")}
\format{
  A data frame with 150 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{IDs of the 150 children}
    \item{\code{o11}}{baseline covariate coded as 0/1: diagnosed with ODD (oppositional defiant disorder) before the first-stage intervention}
    \item{\code{o12}}{baseline covariate: ADHD score at the end of the previous school year (ranging from 0 to 3, larger values for fewer ADHD symptoms)}
    \item{\code{o13}}{baseline covariate coded as 0/1: receiving medication during the previous school year}
    \item{\code{o14}}{baseline covariate coded as 0/1: race - white (coded 1) versus nonwhite (coded 0)}
    \item{\code{a1}}{first-stage intervention coded as -1/1: -1 for low-intensity methamphetamine (MEDS), 1 for low-intensity behavioral modification (BMOD)}
    \item{\code{r}}{first-stage response indicator coded as 0/1}
    \item{\code{o21}}{intermediate outcome: number of months until non-response (maximum: 8 months, NA for responders)}
    \item{\code{o22}}{intermediate outcome coded as 0/1: adherence to the first-stage intervention, 1 for high adherence}
    \item{\code{a2}}{second-stage intervention coded as -1/1: -1 for low-intensity BMOD + MEDS, 1 for high-intensity BMOD}
    \item{\code{y}}{primary outcome (continuous): school performance at the end of the school year (ranging from 1 to 5, higher values reflect better performance)}
  }
}


\references{
Pelham Jr, W. E., & Fabiano, G. A. (2008). Evidence-based psychosocial treatments for attention-deficit/hyperactivity disorder. Journal of Clinical Child & Adolescent Psychology, 37(1), 184-214.
}
\examples{
data(adhd)
attach(adhd)
n = length(a1)
H1 = scale(cbind(o11, o12, o13, o14))
H2 = scale(cbind(H1, a1, H1*a1, r, o22, r*a1, o22*a1))
colnames(H2)[12] = "r*a1"
colnames(H2)[13] = "o22*a1"

fit_ql = ql(H=list(H1, H2), AA=list(a1,a2), RR=list(rep(0, n), y),
            pi=list(rep(0.5, n), rep(0.5,n)), K=2, m=3, lasso=TRUE)

c = 2^c(-15:15)
fit_owl = owl(H=list(H1, H2), AA=list(a1,a2), RR=list(rep(0, n), y),
              n=n, K=2, pi=list(rep(0.5, n), rep(0.5,n)), res.lasso = TRUE,
              loss="hinge", kernel="linear", augment=TRUE, c=c, m=3)


}
\keyword{datasets}
