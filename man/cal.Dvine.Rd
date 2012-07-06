\name{cal.Dvine}
\alias{cal.Dvine}
\title{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines}
\description{Calculating the density of the estimated Dvine at the
  point(s) val.}
\usage{cal.Dvine(obj,val)}
\arguments{
  \item{obj}{object of class 'penDvine', result of 'Dvine'.}
  \item{val}{Values in which the current Dvine should be evaluated.}
}
\details{The current Dvine is evaluated in val and the corresponding
  density values are returned.}
\value{The returing values are the density of the current Dvine at the
  point(s) 'val'.}

\references{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines, Kauermann G. and Schellhase C. (2012), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}