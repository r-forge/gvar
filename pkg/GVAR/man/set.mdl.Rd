\name{set.mdl}
\alias{set.mdl}
\title{Transform VECM to VAR}
\description{
  Transform to the VAR view from a VECM model calculated with \code{est.vecm.mdls} or \code{est.we.mdls}
}
\usage{
set.mdl(mdls, exo = NULL, skip = NULL)
}
\arguments{
  \item{mdls}{an object of class \code{"vecm"}}
  \item{exo}{number of exogenous variables included in the I(1)-model, default is \code{NULL}}
  \item{skip}{skip model parameters that should not be converted to VAR view, default is \code{NULL}}
}

\value{
  A list containing the VECM and VAR views of the model calculated
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

\examples{
##---- Should be DIRECTLY executable !! ----
}

\keyword{ts}

