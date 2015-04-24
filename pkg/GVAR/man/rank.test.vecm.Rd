\name{rank.test.vecm}
\alias{rank.test.vecm}
\alias{rank.test.we}
\title{Test for Cointegrating Rank in Vector Error Correction Models}
\description{
  \code{rank.test.vecm} performs a test for cointegrating rank in VEC models
}
\usage{
rank.test.vecm(Y.ts, etw, p, case, season = NULL, season.start.time = NULL)

rank.test.we(z.ts, etw, p, q = p, n, ex, lex, case)
}
\arguments{
  \item{Y.ts,z.ts}{time series data matrix}
  \item{etw}{estimation time window, a list of the form \code{list(start=,end=,freq=)}}
  \item{p}{lag order of the endogenous variables}
  \item{q}{lag order of the (weakly) exogenous variables, is set to \code{q=p} by default}
  \item{n}{number of endogenous variables}
  \item{ex}{number of exogenous variables}
  \item{lex}{lag order of the strictly exogenous variables}
  \item{case}{intercept and trend options from \code{"I"} to \code{"V"}, where case \code{"I"} is a zero intercept, zero trend model, case \code{"II"} is a restricted intercept, zero trend model, \code{"III"} is a unrestricted intercept, zero trend model, \code{"IV"} is a unrestricted intercept restricted trend model and \code{"V"} is a unrestricted intercept, unrestricted trend model}
  \item{season}{optional seasonal dummies, default value is \code{NULL}}
  \item{season.start.time}{optional dummy start time, default value is \code{NULL}}
}
\details{
  Either use \code{rank.test.vecm} or \code{rank.test.we} for VECM and VECM with exogenous I(1) variables respectively.
}
\value{
  Values for maximum eigenvalue and trace statistics are returned as well as critical values for the limit distribution of these test statistics.
}
\references{Soeren Johansen.
            \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
            Advanced Texts in Econometrics. Oxford University Press, 1995.

            Helmut Luetkepohl, Markus Kraetzig.
            \emph{Applied Time Series Econometrics.}
            Cambridge University Press, 2004.
            
            M. Hashem Pesaran, Yongcheol Shin, and Richard J. Smith.
            Structural analysis of vector error correction models with exogenous I(1) variables.
            \emph{Journal of Econometrics}, 97:293-343, 2000.}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr}
\note{ ~~further notes~~ }
\seealso{\code{\link{est.vecm.mdls}}, \code{\link{est.we.mdls}}}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{ts}
