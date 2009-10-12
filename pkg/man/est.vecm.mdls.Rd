\name{est.vecm.mdls}
\alias{est.vecm.mdls}
\title{Estimation of Vector Error Correction Models}
\description{
  \code{est.vecm.mdls} will calculate VECMs and provide trace and eigenvalue statistics to determine the cointegration rank of the system
}
\usage{
est.vecm.mdls(Y.ts, etw, p, case, r, season = NULL, season.start.time = NULL)
}

\arguments{
  \item{Y.ts}{a time series data matrix containing the endogenous data}
  \item{etw}{estimation time window, a list of the form \code{list(start=,end=,freq=)}}
  \item{p}{lag order of the endogenous variables}
  \item{case}{intercept and trend options from \code{"I"} to \code{"V"}, where case \code{"I"} is a zero intercept, zero trend model, case \code{"II"} is a restricted intercept, zero trend model, \code{"III"} is a unrestricted intercept, zero trend model, \code{"IV"} is a unrestricted intercept restricted trend model and \code{"V"} is a unrestricted intercept, unrestricted trend model}
  \item{r}{cointegrating rank used for estimating the model}
  \item{season}{frequency of optional seasonal dummies, default value is \code{NULL} (no seasonal dummies)}
  \item{season.start.time}{optional dummy start, default value is \code{NULL}}
}
\details{
  Calculating VECMs:
}
\value{
  An object of class \code{"vecm"} containing the following items

  \item{type}{is set to \code{"VECM"} by default}
  \item{dat}{data used for calculating the model}
  \item{freq}{frequency of the time series used}
  \item{n}{number of variables entering the model}
  \item{p}{lag order of the endogenous variables}
  \item{r}{cointegrating rank used for the model estimation}
  \item{T}{length of the time series, given without initial values}
  \item{alpha}{list of model coefficients for the different cointegration ranks, see details}
  \item{beta}{list of model coefficients for the different cointegration ranks, see details}
  \item{Pi}{list of model coefficients for the different cointegration ranks, see details}
  \item{Gamma}{list of model coefficients for the different cointegration ranks, see details}
  \item{case}{intercept and trend specification from \code{"I"} to \code{"V"}, see arguments}
  \item{mu.0}{list of intercepts for the different cointegration ranks, see details}
  \item{mu.1}{list of trend coefficients for the different cointegration ranks, see details}
  \item{Phi}{list of model coefficients for the different cointegration ranks, see details}
  \item{Omega}{list of covariance matrices of the error terms for the different cointegration ranks}
  \item{residuals}{model residuals}
  \item{S}{product moment matrices}
  \item{lambda}{eigenvalues}
  \item{se}{standard errors of coefficients}
  \item{tvalues}{t-distributed test statistic of coefficients}
  \item{pvalues}{p-values for test statistic of coefficients}
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
\note{ ~~further notes~~ }
\seealso{\code{\link{est.we.mdls}}}
\examples{
 data(RDp)

    dt <- 1/tsp(RDp)[3]
     k <- 4
   etw <- list(start=tsp(RDp)[1]+k*dt,end=tsp(RDp)[2],freq=tsp(RDp)[3])

 model <- est.vecm.mdls(Y.ts=RDp,etw=etw,p=k,r=1,case="III",season=4,season.start.time=1)
 summary(model)
}

\keyword{ts}