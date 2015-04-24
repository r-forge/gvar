\name{predict.GVAR}
\alias{predict.GVAR}

\title{Forecasting GVAR objects}
\description{
  forecasting for \code{"GVAR"} objects
}
\usage{
\method{predict}{GVAR}(x, steps, start = NULL, exo.var = NULL, exo.type = "exo")
}

\arguments{
  \item{x}{\code{"GVAR"} object}
  \item{steps}{length of prediction period}
  \item{start}{start prediction at index \code{start}, if \code{NULL} start will be at the end of the sample period}
  \item{exo.var}{either a matrix of dimension exogenous variables x steps if \code{exo.type="exo"} or a vector specifying (subsystem,variable(s))}
  \item{exo.type}{\code{"exo"} for additional exogenous variables, or \code{"from.sub"} to specify in which subsystem the variable is endogenous}
}
\details{
  knows about the future
}
\value{
  Returns a list containing:
  \item{X}{endogenous prediction}
  \item{d}{exogenous prediction}
}
\references{
  GVAR papers
}
\author{
  Rainer Puhr
}

\seealso{
  \code{\link{GVAR}}
}

\keyword{ts}