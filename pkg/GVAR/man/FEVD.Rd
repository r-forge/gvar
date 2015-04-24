\name{FEVD}
\alias{FEVD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forecast Error Variance Decomposition}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FEVD(x, n = 40, shock.var, shock.dir = -1, scal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{shock.var}{
%%     ~~Describe \code{shock.var} here~~
}
  \item{shock.dir}{
%%     ~~Describe \code{shock.dir} here~~
}
  \item{scal}{
%%     ~~Describe \code{scal} here~~
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
function (x, n=40, shock.var, shock.dir=-1,scal=FALSE){

    endoN=sapply(x$we.vecms,function(x) x$n)
    G <- x$G
    U <- x$U
    H.list <- x$H
    H <- H.list[[1]]
    p = x$arguments$p
    q = x$arguments$q
    l <- length(x$subsys)
    if (x$arguments$exo.var) l <- l-1

    if (max(p, q) > 1) {
        for (i in 2:max(p, q)) {
            H <- cbind(H, H.list[[i]])
        }
        I.n <- diag(ncol(H) - ncol(H.list[[max(p, q)]]))
        Zeros <- matrix(0, nrow = ncol(H) - ncol(H.list[[max(p,
            q)]]), ncol = ncol(H.list[[max(p, q)]]))
        H <- rbind(H, cbind(I.n, Zeros))
        G <- rbind(cbind(G, t(Zeros)), cbind(Zeros, I.n))
        U <- rbind(U, matrix(0, nrow(Zeros), ncol(U)))
    }
    G_inv=solve(G)
    Fmat <- G_inv \%*\% H
    UtU= cov(t(U))
#    P <- t(chol(t(UtU)))

    # define shock / either in terms of sd or \%
    s.j <- rep(0, dim(U)[1])
    ind.j <- rep(0, dim(U)[1])
    
    sigma.il <- vector()
      
    for (i in 1:length(shock.dir))
    {    
      if (shock.var[[i]][1] == 1) {
          j <- shock.var[[i]][2]
      }
      else {
          j <- sum(endoN[1:(shock.var[[i]][1] - 1)]) + shock.var[[i]][2]
      }

      # define variables once
      sigma.il[i] <- sqrt(UtU[j, j])

      if(scal){
#        shock=as.numeric(unlist(strsplit(shock.dir,"/sd"))[[1]])*cons*(1/G_inv[j,j])
        shock <- shock.dir[[i]]*sigma.il[i]/(G_inv\%*\%UtU)[j,j]
      }
      else{
        shock <- shock.dir[[i]]
      }
      if(!is.numeric(shock)){
        stop("For the argument shock.dir please submit either a volume shock (e.g. shock.dir='0.003/sd') or a standard deviation shock (e.g. shock.dir='-1').")
      }

      s.j[j] <- shock
      ind.j[j] <- 1
    }

    F.n <- list()
    F.n[[1]] <- diag(nrow(Fmat))
    for (i in 2:(n+1)) {
        F.n[[i]] <- F.n[[i-1]] \%*\% Fmat
    }

    fevd.gi <- matrix(0,length(s.j),n+1)
    fevd.oi <- matrix(0,length(s.j),n+1)
    for (i in which(s.j==0)[1:(sum(endoN)-1)])
    {
      e.i <- rep(0,length(s.j))
      e.i[i] <- 1
      temp <- vector()
      temp2 <- vector()
      temp3 <- vector()
      for (j in 1:(n+1))
      {
        temp[j] <- (t(e.i)\%*\%F.n[[j]]\%*\%UtU\%*\%s.j)^2
        temp2[j] <- t(e.i)\%*\%F.n[[j]]\%*\%UtU\%*\%t(F.n[[j]])\%*\%e.i
#        temp3[j] <- (t(e.i)\%*\%F.n[[j]]\%*\%P\%*\%s.j)^2

        fevd.gi[i,j] <- 1/(t(e.i)\%*\%UtU\%*\%e.i)*sum(temp[1:j])/sum(temp2[1:j])
#        fevd.oi[i,j] <- 1/(t(e.i)\%*\%UtU\%*\%e.i)*sum(temp3[1:j])/sum(temp2[1:j])
      }
    }

    fevd.oi <- NULL
#    psi <- list()
#    psi[[1]] <- t(psi.m[1:endoN[1],])
#    rownames(psi[[1]]) <- 0:n
#    colnames(psi[[1]]) <- colnames(x$Data[[1]])
#    for (i in 2:l)
#    {
#    psi[[i]] <- t(psi.m[(sum(endoN[1:(i-1)])+1):(sum(endoN[1:i])),])
#    rownames(psi[[i]]) <- 0:n
#    colnames(psi[[i]]) <- colnames(x$Data[[i]])
#    }
#    names(psi) <- x$subsys[1:l]

    res <- list(fevd.gi = fevd.gi, fevd.oi = fevd.oi, Fmat = Fmat, G = G, H = H, U = U,shock.dir=shock.dir,sigma.il=sigma.il)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
