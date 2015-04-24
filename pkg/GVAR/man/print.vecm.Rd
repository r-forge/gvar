\name{print.vecm}
\alias{print.vecm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing objects of class vecm}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{vecm}(x, \dots)
}
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  data(RDp)
  p <-  4
  freq <- 4
  etw <- list(start = 1972.25+p*1/freq, end= 1998.75, freq = freq)
  res <- est.vecm.mdls(RDp,etw=etw,p=4,case="III",r=1,season=4,season.start.time=1)
  print(res)
}
\keyword{ts}

