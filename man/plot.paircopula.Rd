\name{plot.paircopula}
\alias{plot.paircopula}
\title{Flexible Pair-Copula Estimation in D-vines with Penalized
  Splines}
\description{Plotting a paircopula of class 'paircopula'.
 }
 \usage{\method{plot}{paircopula}(x,val=NULL,marg=TRUE,plot.view=TRUE,int=FALSE,main.txt=NULL,
                    sub.txt=NULL,contour=FALSE,cuts=20,cex=1,cex.axes=1,
                    xlab=NULL,ylab=NULL,zlab=NULL,xlim=NULL,ylim=NULL,zlim=NULL,
                    show.observ=FALSE,margin.normal=FALSE,...)
}
\arguments{
  \item{x}{object of class 'paircopula', result of function 'paircopula'.}
  \item{val}{Default val = NULL, one can calculate the estimated density/distribution
    for bivariate vector, e.g. val=c(0.5,1).}
  \item{marg}{Default = TRUE, plotting the marginal densities.}
  \item{plot.view}{Default = TRUE, if 'FALSE' no plot is shown, e.g. for
    calculations with val != NULL.}
  \item{int}{Default = FALSE, if TRUE, the integral, i.e. the distribution of
    the copula density is plotted.}
  \item{main.txt}{Default = NULL shows 'K' and the value of lambda.}
  \item{sub.txt}{Default = NULL shows the log-likelihood, the penalized
    log-likelihood and the cAIC-value of the estimation.}
  \item{contour}{If TRUE, a contour plot is shown. Default = FALSE.}
  \item{cuts}{Number of cuts for the contour plots, if contour=TRUE.}
  \item{cex}{Default = 1, determing the size of the main of the plot.}
  \item{cex.axes}{Default = 1, determing the size of the labels at the
    axes.}
  \item{xlab}{Default = NULL and no text is printed at the xlab}
  \item{ylab}{Default = NULL and no text is printed at the ylab}
  \item{zlab}{Default = NULL and 'density' is printed at the zlab for
    int=FALSE and 'distribution' for int=TRUE.}
  \item{xlim}{Default = NULL, changes the range for the values of x in the case of a contour plot.}
  \item{ylim}{Default = NULL, changes the range for the values of y in the case of a contour plot.}
  \item{zlim}{Default = NULL and the range for the values of z are the range of calculated values.}
  \item{show.observ}{Default = FALSE. If TRUE, plotting the original observation
    into a contourplot.}
  \item{margin.normal}{Default = FALSE. If TRUE, the plot is presented with margins following standard normal distribution.}
  \item{...}{further arguments}
}

\value{If 'val' is not NULL, the function returns a matrix with the
  calculated density or distribution values for the set 'val'.
}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann G. and Schellhase C. (2014+), Statistics and Computing (to appear).}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
