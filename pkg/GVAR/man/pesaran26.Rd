\name{pesaran26}
\alias{pesaran26}
\alias{Data}
\alias{weight}
\docType{data}
\title{Pesaran's sample of 26 countries}
\description{
  This dataset contains a set of variables for 26 countries/regions as presented in Pesaran (2007). The variables include real GDP, 
  consumer price index, equity price index, exchange rate, short- and long-term interest rates and an oil price index. Please note 
  that this data is not available for all countries.\cr
  In addition, a matrix containing the average trade weights is provided. 
}
\usage{data(pesaran26)}
\format{
  2 objects can be loaded:
  \tabular{ll}{  
     \code{Data}   \tab a list of length 27, containing the timeseries data of the 26 countries and the exogenous variable.\cr
     \code{weight} \tab a matrix of size 26x26.
  } 
}
\details{
  For further details see Supplement A (Data & Bootstrap) of Pesaran (2007).
}
\source{Stephane Dees, Filippo di Mauro, Hashem Pesaran, and L. Vanessa Smith.
        Exploring the international linkages of the Euro area: A global VAR analysis.
        \emph{Journal of applied Econometrics}, 22(1), 2007.
}
\examples{
data(pesaran26)
# list of countries
names(Data)
# weight matrix
round(weight,3)
}
\keyword{datasets}
