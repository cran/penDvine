\name{my.loop}
\alias{my.loop}

\title{Iterative loop for calculating the optimal coefficients 'v'.}

\description{
  Calculating the optimal coefficients 'v' iteratively, using quadratic programing.
}

\usage{
my.loop(penden.env)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pencopula()}
 }

 \details{
   'my.loop' optimates the log-likelihhod iteratively. Therefore, the
   routine checks a) the relative chance in the optimal lambda and stops the
   iteration, if the relative change of lambda is less than one
   percent. During the calculations of new weights 'v' in the routine
   'new.weights', most of the values are called '.temp'. This add on
   underlines the temporarily values. Alternatively b) for fixed lambda, 'my.loop' checks the relative
   change in the weights. If the change of a) the optimal lambda or b)
   of the basis coefficients 'v' are greater than one percent, the
   the real values are overwritten with the '.temp' values.
}

\value{
  \item{liste}{The results of each iteration are written in a matrix called 'liste',
   saved in the environment. 'liste' contains the penalized
   log-likelihood, the log-likelihood, 'lambda' and the weights 'v'.}
}

\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{math}


