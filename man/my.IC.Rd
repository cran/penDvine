\name{my.IC}
\alias{my.IC}

\title{Calculating the AIC-, cAIC- and BIC-value}

\description{
Calculating the AIC-, cAIC- and BIC- value of the paircopula density estimation. Therefore, we add the unpenalized log likelihood of the estimation and the degree of freedom. 	
}
\usage{
my.IC(penden.env,temp=FALSE)
}

\arguments{
  \item{penden.env}{Containing all information, environment of
    paircopula()}
  \item{temp}{Default=FALSE, if TRUE temporary values of AIC, cAIC and
  BIC are calculated.}
  
}

\details{
AIC is calculated as
\eqn{AIC(\lambda)= - 2*l({\bf u},\hat{\bf{v}}) + 2*df(\lambda)}

cAIC is calculated as
\eqn{AIC(\lambda)= - 2*l({\bf u},\hat{\bf{v}}) + 2*df(\lambda)+(2*df*(df+1))/(n-df-1)}

BIC is calculated as
\eqn{BIC(\lambda)= 2*l({\bf u},\hat{\bf{v}}) + 2*df(\lambda)*log(n)}
}

\value{
  \item{AIC}{sum of twice the negative non-penalized log likelihood and
    mytrace}
  \item{cAIC}{corrected AIC.}
  \item{trace}{calculated mytrace as the sum of the diagonal matrix
    df, which results as the product of the inverse of the penalized
    second order derivative of the log likelihood with the non-penalized
    second order derivative of the log likelihood}
  \item{BIC}{sum of twice the non-penalized log likelihood and log(n)}
  All values are saved in the environment.
}
\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{math}

