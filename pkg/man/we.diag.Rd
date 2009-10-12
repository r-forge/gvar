\name{we.diag}
\alias{we.diag}
\title{Exogeneity diagnostics}
\description{
  This function tests variables for weak exogeneity.
}
\usage{
we.diag(Z, n, p, q = p, case)
}

\arguments{
  \item{Z}{a time series data matrix}
  \item{n}{number of endogenous variables, the remaining variables will be tested for weak exogeneity}
  \item{p}{lag order of the endogenous variables}
  \item{q}{lag order of the weakly exogenous variables}
  \item{case}{intercept and trend options from \code{"I"} to \code{"V"}, where case \code{"I"} is a zero intercept, zero trend model, case \code{"II"} is a restricted intercept, zero trend model, \code{"III"} is a unrestricted intercept, zero trend model, \code{"IV"} is a unrestricted intercept restricted trend model and \code{"V"} is a unrestricted intercept, unrestricted trend model}
}
\details{
  In order to be weakly exogenous variables X must be (a) integrated of order 1 (X~I(1)); it follows (b) that X is not cointegrated on its own, and (c) that the differenced process does not depend on the lagged Z.
  
  The variables tested will be the variables from the last (N-n) columns of the matrix \code{Z}, where N is the total number of columns.
}
\value{
  Results for test (a) to (c).
}
\references{Stefan Zeugner.
            \emph{Implementing Pesaran-Shin-Smith.}
            First year paper, Institute for Advanced Studies, Vienna, 2006.

            Soeren Johansen.
            \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
            Advanced Texts in Econometrics. Oxford University Press, 1995.

            M. Hashem Pesaran, Yongcheol Shin, and Richard J. Smith.
            Structural analysis of vector error correction models with exogenous I(1) variables.
            \emph{Journal of Econometrics}, 97:293-343, 2000.}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr}

\seealso{\link{rank.test.vecm}}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{ts}

