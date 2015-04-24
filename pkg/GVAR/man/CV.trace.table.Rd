\name{CV.trace.table}
\alias{CV.maxeigen.table}
\docType{data}
\title{Estimates of critical values of the limiting distributions of the cointegrating rank statistics in a VECM with weakly exogenous I(1) variables}
\description{
  These data sets provide estimated critical values for the limiting distributions of the cointegrating rank trace (\code{CV.trace.table}) and maximum eigenvalue (\code{CV.maxeigen.table}) statistics for all five intercept/trend scenarios (case I to V) as describe in PSS 2000.
}
\usage{data(cv.tables.pesaran)}
\format{
  The format is:
Lists of 10
 $ case I 5%
 $ case I 10%
 $ case II 5%
 $ case II 10%
 $ case III 5%
 $ case III 10%
 $ case IV 5%
 $ case IV 10%
 $ case V 5%
 $ case V 10%
}
\details{
  details to follow
}
\source{
  M. Hashem Pesaran, Yongcheol Shin, and Richard J. Smith.
  Structural analysis of vector error correction models with exogenous I(1) variables.
  \emph{Journal of Econometrics}, 97:293-343, 2000.
  
  James G. MacKinnon, Alfred A. Haug and Leo Michelis.
  Numerical Distribution Functions of Likelihood Ratio Tests for Cointegration.
  \emph{Journal of Applied Econometrics}, 5(14):563-577, 1999.
}

\examples{
data(CV.trace.table)
data(CV.maxeigen.table)

CV.trace.table[[1]][12,1]         ## returns the 0.05 trace CV for case I in a ... model
CV.maxeigen.table[[4]][12,1]      ## returns the 0.1 maximum eigenvalue CV for case II in a ... model
}
\keyword{datasets}
