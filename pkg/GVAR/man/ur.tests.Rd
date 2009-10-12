\name{ur.tests}
\alias{ur.tests}
\title{Unit Root Tests}
\description{
  \code{ur.tests} computes the three different unit root test statistics, Augmented Dickey-Fuller Test, Phillips-Perron Unit Root Test and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) Test, from \pkg{tseries} for a data matrix.
}
\usage{
ur.tests(X)
}

\arguments{
  \item{X}{a (time series) data matrix}
}

\value{
  Returns a matrix giving the alternative hypothesis of the tests and the p-values for each variable.
}
\references{see package \pkg{tseries}}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr}

\examples{
X <- ts(cbind(rnorm(100),rnorm(100),rnorm(100)),start=1980,freq=4)
ur.tests(X)

data(RDp)
ur.tests(RDp)
}

\keyword{ts}

