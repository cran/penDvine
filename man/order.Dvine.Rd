\name{order.Dvine}
\alias{order.Dvine}
\title{Ordering the first level of the Dvine.}

\description{
Ordering the first level of the Dvine, depending on the pairwise cAIC-values.
}

\usage{
order.Dvine(help.env)
}

\arguments{
  \item{help.env}{Containing all information, environment of
    Dvine()}
}

\details{
Beginning in the top tree level of a p-dimensional D-vine, we calculate all
(p,2) 'p over 2' marginal pairwise copulas fitted by
penalized splines. For each pair (i,j) this gives the fitted maximized
likelihood value. We order the variable pairs, subject
to their increasing estimated pairwise cAIC and start with the pair of
covariates with lowest estimated cAIC. We now select the pairs of
variables such that the resulting selection gives a tree. The problem of finding this
selection is equivalent to solve a traveler salesman problem by interpreting the cAIC as distance measure between
two variables. Once this problem is
solved, the specification of the first tree level completely defines the
D-vine
}

\value{
  \item{order}{Calculated order of the first level of the Dvine, saved
  in the environment.}
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized Splines, Kauermann G. and Schellhase C. (2012), working paper}
\keyword{nonparametric}
