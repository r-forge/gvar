\name{est.we.mdls}
\alias{est.we.mdls}
\title{Estimation of Weakly Exogenous VEC Models}
\description{
  \code{est.we.mdls} will calculate weakly exogneous VECMs and provide trace and eigenvalue statistics to determine the cointegration rank of the system
}
\usage{
est.we.mdls(z.ts, etw, p, q = p, n, ex=0, lex=NULL, case, r)
}
\arguments{
  \item{z.ts}{a time series data matrix containing endogenous and weakly exogenous variables}
  \item{etw}{estimation time window, a list of the form \code{list(start=,end=,freq=)}}
  \item{p}{lag order of endogenous variables}
  \item{q}{lag order of (weakly) exogenous variables}
  \item{n}{number of endogenous variables}
  \item{ex}{number of exogenous variables, set to \code{0} by default}
  \item{lex}{lag order of exogenous variables, set to \code{NULL} by default}
  \item{case}{intercept and trend options from \code{"I"} to \code{"V"}, where case \code{"I"} is a zero intercept, zero trend model, case \code{"II"} is a restricted intercept, zero trend model, \code{"III"} is a unrestricted intercept, zero trend model, \code{"IV"} is a unrestricted intercept restricted trend model and \code{"V"} is a unrestricted intercept, unrestricted trend model}
  \item{r}{cointegrating rank used for model estimation}
}
\details{Calculating weakly exogneous VECMs:
}
\value{
  An object of class \code{"vecm"} containing the following items

  \item{type}{is set to \code{"weakly exogenous VECM"} by default}
  \item{dat}{data used for calculating the model}
  \item{freq}{frequency of the time series used}
  \item{m}{number of total variables entering the model}
  \item{n}{number of endogenous variables entering the model}
  \item{p}{lags for the endogenous variables}
  \item{q}{lags for the (weakly) exogenous variables}
  \item{r}{cointegrating rank used for model estimation}
  \item{T}{length of the time series, given without initial values}
  \item{alpha}{list of model coefficients for the different cointegration ranks, see details}
  \item{beta}{list of model coefficients for the different cointegration ranks, see details}
  \item{Pi.y}{list of model coefficients for the different cointegration ranks, see details}
  \item{Phi}{list of model coefficients for the different cointegration ranks, see details}
  \item{Psi}{list of model coefficients for the different cointegration ranks, see details}
  \item{Lambda}{list of model coefficients for the different cointegration ranks, see details}
  \item{case}{intercept and trend specification from \code{"I"} to \code{"V"}, see arguments}
  \item{c.0}{list of intercepts for the different cointegration ranks, see details}
  \item{c.1}{list of trend coefficients for the different cointegration ranks, see details}
  \item{Omega.uu}{list of covariance matrices of the error terms for the different cointegration ranks}
  \item{S}{product moment matrices}
  \item{lambda}{eigenvalues}
  \item{residuals}{model residuals}
}
\references{Soeren Johansen.
            \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
            Advanced Texts in Econometrics. Oxford University Press, 1995.

            James G. McKinnon et al.
            \emph{Numerical Distribution Functions of Likelihood Ratio Test for Cointegration.}
            

            M. Hashem Pesaran, Yongcheol Shin, and Richard J. Smith.
            Structural analysis of vector error correction models with exogenous I(1) variables.
            \emph{Journal of Econometrics}, 97:293-343, 2000.}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr}
\note{ ~~further notes~~

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{\code{\link{est.vecm.mdls}}}
\examples{
##---- Should be DIRECTLY executable !! ----
}

\keyword{ts}