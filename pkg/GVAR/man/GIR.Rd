\name{GIR}
\alias{GIR}
\title{Generalized Impulse Response}
\description{Perform Generalized Impulse Response analysis.}
\usage{
GIR(x, n, shock.var)
}
\arguments{
  \item{x}{GVAR objekt}
  \item{n}{time interval from T to T+n}
  \item{shock.var}{vector of length 2 containing the index of the country and the variable to be shocked}
}
\details{to be added}
\value{a list containing:
  \item{psi}{impulse responses}
  \item{Fmat}{invG*H}
  \item{G}{multiplier matrix for the current variables}
  \item{H}{multiplier matrix for the lagged variables}
  \item{U}{residuals}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Rainer Puhr, Martin Feldkircher}
\seealso{\code{\link{GVAR}}}
\examples{}
\keyword{ts}