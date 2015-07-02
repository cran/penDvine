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
Version: \tab 0.2.4\cr
Date: \tab 2015-07-02\cr
License: GPL (>= 2)
LazyLoad: \tab yes\cr
}
The packages contributes the function 'Dvine()' for estimating densities
and distributions of Dvines using penalized splines techniques. 
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014), Statistics and Computing 24(6): 1081-1100).}
\examples{
#This examples describes the estimation of a D-vine to winddata,
#available in this package. After the margins are prepared, we estimate a
#D-vine using B-splines with 9 (K+1) marginal knots and penalizing second (m=2)
#order differences (pen=1) of the basis coefficients.

data(Winddata)
ex.data<-Winddata[c(1:100),c(1:4)] #examplary subdata for fast calculation
wind.example<-Dvine(ex.data,K=7,pen=1,lambda=100,order.Dvine=FALSE,base="B-spline",m=2,cores=2)
  }

	   
	   
