\name{print.vecm}
\alias{print.vecm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.vecm(x, ...)
}
%- maybe also 'usage' for other objects documented here.
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

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, ...) 
{
    if (x[["type"]] == "pure VECM") {
        roundto <- 5
        alpha <- round(x$alpha, roundto)
        beta <- round(x$beta, roundto)
        rownames(beta) <- colnames(x$dat)
        colnames(beta) <- paste(1:x$r, ".:   ", sep = "")
        Gamma <- x$Gamma
        const <- cbind(x$mu0, x$mu1)
        nam <- vector()
        if (!is.null(x$mu0)) {
            nam <- c(nam, "         Const")
        }
        if (!is.null(x$mu1)) {
            nam <- c(nam, "Trend")
        }
        if (!is.null(const)) {
            const <- round(const, roundto)
            colnames(const) <- nam
        }
        if (!is.null(x$season)) {
            Phi <- round(x$Phi, roundto)
            nam <- "            s1"
            if (x$season > 1) {
                nam <- c(nam, paste("s", 2:(x$season - 1), sep = ""))
            }
            colnames(Phi) <- nam
        }
        cat("\n")
        cat("Coefficients:\n")
        cat("Beta:\n")
        print(t(beta))
        cat("\n")
        cat("Alpha:\n")
        print(alpha)
        cat("\n")
        if (x$p > 1) {
            Gammer <- NULL
            cat("Gamma:\n")
            for (i in 1:length(x$Gamma)) {
                colnames(Gamma[[i]])[1] <- paste("        ", 
                  colnames(Gamma[[i]])[1])
                Gammer <- cbind(Gammer, round(Gamma[[i]], roundto))
            }
            print(Gammer)
            cat("\n")
        }
        if (!is.null(x$mu0) || !is.null(x$mu1)) {
            cat("Intercept (and Trend) in VAR:\n")
            print(const)
            cat("\n")
        }
        if (!is.null(x$season)) {
            cat("Seasonality:\n")
            print(Phi)
            cat("\n")
        }
        cat("\n")
    }
    else if (x[["type"]] == "weakly exogenous VECM") {
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
