\name{new.weights}
\alias{new.weights}
\title{Calculating new weights v.}

\description{
Calculating new weights v using quadratic programing.
}

\usage{
new.weights(penden.env, lambda.temp=NULL)
}

\arguments{
  \item{penden.env}{Containing all information, environment of
    paircopula()}
  \item{lambda.temp}{Default=NULL, if optimal lambda is calculated, the
  lambda are saved temporarily, also the resulted coefficients are saved
  temporarily until some convergences.
  }
}

\details{
  The new weights are calculated solving a quadratic program. Therefore,
  the derivates of first and second order are needed, 'Derv1.pen' and
  'Derv2.pen'. Moreover, we have to fulfill the side conditions
  v>=0, sum(v)=1 and that the marginal densities are uniform. All side
  conditions are saved as 'AA.help' in the environment. 

  If the quadratic program does not find a new feasible solution, the whole
  program terminates. For solving the quadratic program, we use the
  function 'solve.QP' from the  R-package 'quadprog'.
}

\value{
  \item{ck.val.temp}{Calculated new values for the weights 'v'. The add
    on 'temp' means, that there is a check in the next step if the
    weights 'v' have been converted (in the case of fixed lambda). If converted, the new values
    'ck.val.temp' are unnoted. If not converted, 'ck.val.temp' become
    the ordinary 'ck.val' for the next iteration. This check is done in
    my.loop.
    If the optimal value of lambda is calculated, the coefficients 'ck.val.temp' become
    the ordinary 'ck.val' for the next iteration if lambda is converted.
    t
  }
  'ck.val.temp' is saved in the environment.
}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014+), Statistics and Computing (to appear).}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{nonparametric}
