\name{summary.vecm}
\alias{summary.vecm}
\title{Summarizing objects of class vecm}
\description{\code{summary.vecm} plots summary statistics for objects of class \code{vecm} 
}
\usage{
\method{summary}{vecm}(object, ...)
}

\arguments{
  \item{object}{object of class \code{vecm}}
  \item{\dots}{}
}
\details{
  \code{summary.vecm} will provide an overview of the coefficients, standard errors etc. of an \code{vecm}-object 
}
\value{

}
\references{Soeren Johansen.
            \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
            Advanced Texts in Econometrics. Oxford University Press, 1995.

            Helmut Luetkepohl, Markus Kraetzig.
            \emph{Applied Time Series Econometrics.}
            Cambridge University Press, 2004.
            
            Helmut Luetkepohl.
            \emph{New Introduction to Multiple Time Series Analysis.}
            Springer, 2005.}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr, Michael Sigmund}            
\note{
}
\seealso{\code{\link{est.vecm.mdls}}, \code{\link{print.vecm}}}
\examples{
  data(RDp)
  p <-  4
  freq <- 4
  etw <- list(start = 1972.25+p*1/freq, end= 1998.75, freq = freq)
  res <- est.vecm.mdls(RDp,etw=etw,p=4,case="III",r=1,season=4,season.start.time=1)
  summary(res)
}

\keyword{ts}

