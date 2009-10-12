\name{RDp}
\alias{RDp}
\docType{data}
\title{German long term interest and inflation rate}
\description{
  This data set gives the quarterly German nominal long term interest rate (source: Monatsberichte der Deutschen Bundesbank, quarterly values are values of last month of quarter) and the quarterly German inflation rate (source: Deutsches Institut für Wirtschaftsforschung, Volkswirtschaftliche Gesamtrechnung). The sample period is the second quarter of 1972 to the fourth quarter of 1998.
}
\usage{data(RDp)}
\format{
  A time series matrix containing 2 columns of 107 observations each.
}
\details{
  The data before 1990Q3 is West Germany data only, all of Germany afterwards.
}
\source{Helmut Luetkepohl.
        \emph{New Introduction to Multiple Time Series Analysis.}
        Springer, 2005.}
\examples{
data(RDp)
plot(RDp)
}
\keyword{datasets}
