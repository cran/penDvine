\name{bernstein}
\alias{bernstein}
\alias{bernstein2}
\title{Flexible Pair-Copula Estimation in D-vines with Penalized Splines}
\description{Calculation of Bernstein Polynomials, following the formula
(choose(n, v) * x^v * (1 - x)^(n - v)) * (n + 1) .}
\usage{bernstein(v,x,n)}
\arguments{
  \item{v}{n choose v.}
  \item{x}{Argument between 0 and 1.}
  \item{n}{n choose v.}
}
\note{bernstein and bernstein2 are identical, only the first argument of
the functions differs for different applications of 'apply'.}
\value{Returning a matrix, containing Bernstein Polynomials.}
\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
