\name{Winddata}
\alias{Winddata}
\alias{Winddata-original}
\docType{data}
\title{Maximal daily windspeed in khm.}
\description{Maximal daily windspeed in khm at Frankfurt, Berlin, Bremen,
  Munich and Cuxhaven, observed in the time from 01/01/2000 until 31/12/2011.}
\usage{data(Winddata)}

\format{
A data frame with 507 observations of the following 2 variables.
  \describe{
	\item{\code{Date}}{Date}
	\item{\code{Frankfurt}}{Observations in Frankfurt}
	\item{\code{Berlin}}{Observations in Berlin}
	\item{\code{Bremen}}{Observations in Bremen}
	\item{\code{Munich}}{Observations in Munich}
	\item{\code{Cuxhaven}}{Observations in Cuxhaven}
  }
}

\note{The data is available at the internet page of the 'Deutsche
  Wetterdienst' DWD, www.dwd.de.
  'Winddata-original' contains the original data, 'Winddata' contains
  marginal fitted t-distributions of the data.
}

\keyword{datasets}