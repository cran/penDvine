\name{Derv1}
\alias{Derv1}
\title{Calculating the first derivative of the paircopula likelihood function w.r.t. parameter b}
\description{
 Calculating the first derivative of the paircopula likelihood function
 w.r.t. parameter v. 
}
\usage{
Derv1(penden.env,temp=FALSE,lambda=NULL)
}
\arguments{
	\item{penden.env}{Containing all information, environment of
	  paircopula().}
	\item{temp}{Default=FALSE,if TRUE temporary calculations of
	  optimal parameters are done.}
	\item{lambda}{Default=NULL, i.e. the saved smoothing parameter lambda in the
	  environment is used. Alternatively, temporary values of lambda
	  are used for optimization of lambda.}
      }

\details{
The calculation of the first derivative of the paircopula likelihood function w.r.t. b equals
\deqn{s(v,\lambda)= {\partial l(v,\lambda)}/{\partial v}= \sum_{i=1}^n
  \Phi(u_i)/c(u_i,v) - P(\lambda)v}{%
  \eqn{s(v,lambda)=}}
  with \deqn{P(\lambda)}{%
    \eqn{P(lambda)}}
  is the penalty matrix, saved in the environment.
}

\value{
  \item{Derv1.pen}{first order derivation of the penalized likelihood
    function w.r.t. parameter v.}
  Derv1.pen is saved in the environment.
}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014+), Statistics and Computing (to appear).}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{math}
