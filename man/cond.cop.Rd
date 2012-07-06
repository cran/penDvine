\name{cond.cop}
\alias{cond.cop}
\title{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines}
\description{Calculation of the conditional paircopulas.}
\usage{cond.cop(data,coef,K,diff="u2",Index.basis.D,base,q=NULL)}
\arguments{
  \item{data}{Considered bivariate data, later used as "u1" and "u2".}
  \item{coef}{Considered coefficients of the splines.}
  \item{K}{Number of splines.}
  \item{diff}{Default="u2", alternatively diff="u1". Determining in which
    direction the expression is differentiated.}
  \item{Index.basis.D}{Vector of indices built in the program before.}
  \item{base}{"B-spline" or "Bernstein".}
  \item{q}{Order of the B-spline, default="NULL"}
}

\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
