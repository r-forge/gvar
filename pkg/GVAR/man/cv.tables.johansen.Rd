\name{cv.tables.johansen}
\alias{cv.tables.johansen}
\docType{data}
\title{Estimates of critical values of the limiting distributions of the cointegrating rank statistics in a VECM}
\description{
  This data set provides estimated critical values for the limiting distributions of the cointegrating rank trace (\code{CV.trace.table}) and maximum eigenvalue (\code{CV.maxeigen.table}) statistics for all five intercept/trend scenarios (case I to V) as describe in Joh 1995.
}
\usage{data(cv.tables.johansen)}
\format{
List of 5 tables. Each table gives the 90\%-, 95\%- and 99\%-quantiles for the cointegration rank test for \code{p} variables minus rank \code{r}.
}
\details{
  \code{table15.1} gives the tables for the zero intercept, zero trend model, \code{table15.2} for the restricted intercept, zero trend model, \code{table15.3} for the unrestricted intercept, zero trend model, \code{table15.4} for the unrestricted intercept, restricted trend model, and \code{table15.5} for the unrestricted intercept unrestricted trend model.
}
\source{Soeren Johansen.
        \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
        Advanced Texts in Econometrics. Oxford University Press, 1995.
}

\examples{
data(cv.tables.johansen)
cv.tables.johansen$table15.1
}
\keyword{datasets}
