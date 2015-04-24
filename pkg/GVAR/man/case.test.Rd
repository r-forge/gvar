\name{case.test}
\alias{case.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test case
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
case.test(z.ts, etw, p, q = NULL, n = NULL, r, type = "vecm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z.ts}{
%%     ~~Describe \code{z.ts} here~~
}
  \item{etw}{
%%     ~~Describe \code{etw} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(z.ts,etw,p,q=NULL,n=NULL,r,type="vecm")
   #  z.ts ... time series
   #   etw ... estimation time window as specified in estim.we.mdls.r
   #     p ... lag order
   #     n ... number of endogenous variables
   #     r ... cointegrating relations
   #  type ... model, either "vecm" or "we.vecm"
{
  if (is.null(n) && type=="we.vecm") stop("\'n\' missing",call.=FALSE)
  
  mdls <- list()
  test.stat <- vector()
  cv <- vector()
  p.value <- vector()
  
  T <- dim(z.ts)[1]-p
  m <- dim(z.ts)[2]
  cases <- c("I","II","III","IV","V")
  
  if (type=="we.vecm")
  {
  
     #     Johansen      ~      Pesaran
   
     # case<- "H_2(r)"   ~  case "I"    # \mu_t = 0
     # case<- "H_1^*(r)" ~  case "II"   # \mu_t = \alpha \rho_0
     # case<- "H_1(r)"   ~  case "III"  # \mu_t = \mu_0
     # case<- "H^*(r)"   ~  case "IV"   # \mu_t = \mu_0+\alpha \rho_1 t # \alpha_{\bot}'\mu_1=0
     # case<- "H(r)"     ~  case "V"    # \mu_t = \mu_0+\mu_1 t
      
    for (i in 1:5){mdls[[i]] <- est.we.mdls(z.ts=z.ts,etw=etw,p=p,q=q,r=r,n=n,case=cases[i])}

  } else {
    for (i in 1:5){mdls[[i]] <- est.vecm.mdls(Y.ts=z.ts,etw=etw,p=p,r=r,case=cases[i])}
  }
  
    # test statistics as in Johansen (1995) p. 161 f.
  
  for (i in c(1,3))
  {
    test.stat[i] <- T*sum(log(1-mdls[[5-i+1]]$lambda[(r+1):m])-log(1-mdls[[4-i+1]]$lambda[(r+1):m]))
    cv[i] <- qchisq(0.95,df=(m-r))
    p.value[i] <- 1-pchisq(test.stat[i],df=(m-r))
  }
  for (i in c(2,4))
  {
    test.stat[i] <- T*sum(log(1-mdls[[3-i+2]]$lambda[1:r])-log(1-mdls[[4-i+2]]$lambda[1:r]))
    cv[i] <- qchisq(0.95,df=r)
    p.value[i] <- 1-pchisq(test.stat[i],df=r)
  }

    # prepare & print output
      
  res <- cbind(round(test.stat,4),round(cv,4),round(p.value,4))
  rownames(res) <- c("c.1 restricted vs. c.1 unrestricted","       c.1 = 0 vs. c.1 restricted",
                     "c.0 restricted vs. c.0 unrestricted","       c.0 = 0 vs. c.0 restricted")
  colnames(res) <- c("     Statistic","   Crit.value","    p-value")
  
  cat("\nCase selection results:\n")
  print(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
