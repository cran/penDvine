\name{pen.log.like}
\alias{pen.log.like}
\title{Calculating the log likelihood}

\description{
Calculating the considered log likelihood.
}

\usage{
pen.log.like(penden.env, cal=FALSE, temp=FALSE)
}

\arguments{
  \item{penden.env}{Containing all information, environment of paircopula()}
  \item{cal}{if TRUE, the final weights of one iteration are used for
    the calculation of the penalized log likelihood.}
  \item{temp}{if TRUE, the iteration for optimal weights is still in
    progress and the temporary weights are used for calculation.}
}
\details{
The calculation depends on the estimated weights v, the penalized
splines Phi and the penalty paramters lambda.\cr
\deqn{l(v,\lambda)=\sum_{i=1}^{n} \left[ \log \{\sum_{i=1}^n
  \boldsymbol\Phi(u_i)\} v\right]- \frac 12 v^T \boldsymbol{P}(\lambda) b}{%
  \eqn{l(v,lambda)=sum(log(\Phi(u_i)v))-0.5*b^T \tilde{P} v}}

The needed values are saved in the environment.
}

\value{
  \item{pen.log.like}{Penalized log likelihood of the paircopula density.}
  \item{log.like}{Log-Likelihood of the paircopula density.}
The values are saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}

\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014+), Statistics and Computing (to appear).}
