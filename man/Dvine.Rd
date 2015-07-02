\name{Dvine}
\alias{Dvine}
\title{Flexible Pair-Copula Estimation in D-vines with Penalized Splines}
\description{Calculating D-vines with penalized B-splines or penalized
  Bernstein polynomials}
\usage{Dvine(data,K=8,lambda=100,order.Dvine=TRUE,pen=1,base="Bernstein",m=2,cores=NULL,q=2)}
     
\arguments{
  \item{data}{'data' contains the data. 'data' has to be a matrix or a
    data.frame with two columns.}
  \item{K}{K is the degree of the Bernstein polynomials. In the case of
    linear B-spline basis functions, K+1 nodes are used for the basis functions.}
  \item{lambda}{Starting value for lambda, default is lambda=100.}
  \item{order.Dvine}{Indicating if the first level of the Dvine is
    ordered, default (order.Dvine=TRUE).}
  \item{pen}{'pen' indicates the used penalty. 'pen=1' for the
    difference penalty of m-th order. 'pen=2' is only implemented for
    Bernstein polynomials, it is the penalty based
    on the integrated squared second order derivatives of the Bernstein polynomials.}
  \item{base}{Type of basis function, default is
    "Bernstein". An alternative is base="B-spline".}
  \item{m}{Indicating the order of differences to be penalised. Default is "m=2".}
  \item{cores}{Default=NULL, the number of cpu cores used for parallel computing can be specified.}
  \item{q}{Order of B-splines. Default is q=2, NULL if Bernstein polynomials are used.}
}
\details{The calculation of the Dvine is done stepwise. If the option
  'order.Dvine' is selected, the order of the first level of the Dvine
  is specifed. From the second level, each paircopula is calculated
  (parallel or not) until the highest level.
  The specifications in 'Dvine' are done for every paircopula in the
  Dvine. There is no option to change parameters for some paircopulas.  
}
\value{Returning a list containing
  \item{Dvine}{an object of class 'Dvine'}
  \item{log.like}{the estimated log-likelihood}
  \item{AIC}{AIC value}
  \item{cAIC}{corrected AIC value}
  \item{K}{Number of K}
  \item{order}{the used order of the first level}
  \item{S}{Sequence seq(1:(dim(data)[2]))}
  \item{N}{Number of observations, that is dim(data)[1]}
  \item{base}{Used basis function}
  }

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014+), Statistics and Computing (to appear).}
