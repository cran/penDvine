\name{penDvine-package}
\alias{penDvine-package}
\docType{package}
\title{
The package 'penDvine' offers routines for estimating densities and
copula distribution of D-vines.
}
\description{
The package 'penDvine' offers routines for estimating densities and
distribution of D-vines. For details see the description of the function Dvine().
}
\details{
\tabular{ll}{
Package: \tab penDvine\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2012-07-19\cr
License: GPL (>= 2)
LazyLoad: \tab yes\cr
}
The packages contributes the function 'Dvine()' for estimating densities
and distributions of Dvines using penalized splines techniques. 
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\examples{
  #This examples describes the estimation of a D-vine to winddata,
#available in this package. After the margins are prepared, we estimate a
#D-vine using B-splines with 9 (K+1) marginal knots and penalizing second (m=2)
#order differences (pen=0) of the basis coefficients.

  data(Winddata)
  print(dim(Winddata))
  bspl.diff2 <-  Dvine(Winddata,K=8,pen=1,lambda=1000,order.Dvine=FALSE,base="B-spline",m=2)
  }

	   
	   