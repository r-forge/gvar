\name{GVAR}
\alias{GVAR}
\title{Global Vector Auto-Regressive Modelling}
\description{
  \code{GVAR} computes VECMs for all regions and stacks the models to a Global Vector Autoregressive Model
}
\usage{
GVAR(data, tw = NULL, p, q = p, r = NULL, weight, case, exo.var = FALSE,
     d = NULL, endo = NULL, ord = NULL, we = NULL, method = "max.eigen")
}

\arguments{
  \item{data}{timeseries data as list (each entry is a matrix of a subsystem of variables, if \code{exo.var=TRUE} the last entry are exogeneous variables)}
  \item{tw}{time window, vector of start and end point, if \code{NULL} the maximum time interval will be used}
  \item{p}{scalar/vector of endogenous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{q}{scalar/vector of weakly exogeneous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{r}{scalar/vector of cointegrating relations, if a scalar is provided the same cointegration rank is used for all subsystems, if set to \code{NULL} the cointegration rank will be determind by \code{method}}
  \item{weight}{weight matrix, the diagonal elements need to be zero}
  \item{case}{single value/vector of cases (\code{"I"} to \code{"V")}, where case \code{"I"} is a zero intercept, zero trend model, case \code{"II"} is a restricted intercept, zero trend model, \code{"III"} is a unrestricted intercept, zero trend model, \code{"IV"} is a unrestricted intercept restricted trend model and \code{"V"} is a unrestricted intercept, unrestricted trend model; if a single value is provided the same structure is used for all subsystems}
  \item{exo.var}{if \code{TRUE} strictly exogeneous variables are included in the model}
  \item{d}{list showing which strictly exogeneous variables enter the subsystem equations, if \code{NULL} all variables in the last entry of \code{data} will be used for all subsystems}
  \item{lex}{scalar/vector of exogeneous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{endo}{list of endogenous variables used in each subsystem, if \code{NULL} all variables in \code{data} will be used}
  \item{ord}{vector used if variables in the different subsystem don't appear in the same order, order of each subsystem is concatenated to one vector, if \code{NULL} the variables in \code{data} are assumed to be ordered identically in all subsystems}
  \item{we}{list with numbers of weakly exogeneous variables included in each VECM, corresponds to numbers in \code{ord}, if \code{NULL} every variable appearing in all subsystems will be used}
  \item{method}{select cointegrating rank by max. eigenvalue (\code{"max.eigen"}) or trace statistic (\code{"trace"})}
}
\details{
  The function computes a VECM for every subsystem before stacking the results to a GVAR model.
  
  Specification of input here.
}
\value{
  An object of class \code{"GVAR"} containing the following items:
  
  \item{subsys}{subsystem names}
  \item{G}{multiplier matrix for the current variables}
  \item{H}{multiplier matrix for the lagged variables}
  \item{Upsilon.0}{multiplier matrix for the current strictly exogenous variables}
  \item{Upsilon}{multiplier matrix for the lagged strictly exogenous variables}
  \item{c.0}{multiplier matrix for the intercept}
  \item{c.1}{multiplier matrix for the trend}
  \item{we.vecms}{VECMs of the subsystems}
  \item{m}{total variables in the subsystems}
  \item{n}{endogenous variables in the subsystems}
  \item{p}{endogenous lags in the subsystems}
  \item{q}{(weakly) exogenous lags in the subsystems}
  \item{r}{cointegration rank in the subsystems}
  \item{weight}{weight matrix used to calculate the weakly exogenous variables}
  \item{U}{residuals of the GVAR}
  \item{U.cov}{residual covariance matrix}
}
\references{Stephane Dees, Filippo di Mauro, Hashem Pesaran, and L. Vanessa Smith.
            Exploring the international linkages of the Euro area: A global VAR analysis.
            \emph{Journal of applied Econometrics}, 22(1), 2007.

            Soeren Johansen.
            \emph{Likelihood-Based Inference in Cointegrated Vector Auto-Regressive Models.}
            Advanced Texts in Econometrics. Oxford University Press, 1995.

            M. Hashem Pesaran, Yongcheol Shin, and Richard J. Smith.
            Structural analysis of vector error correction models with exogenous I(1) variables.
            \emph{Journal of Econometrics}, 97:293-343, 2000.}
\author{Martin Summer, Klaus Rheinberger, Rainer Puhr}

\seealso{\code{\link{est.we.mdls}}}

\examples{}

\keyword{ts}

