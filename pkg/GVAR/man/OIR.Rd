\name{OIR}
\alias{OIR}
\title{
Orthogonalized Impulse Response
}
\description{
Orthogonalized Impulse Response for GVAR objects
}
\usage{
OIR(x, n = 40, shock.var, shock.dir = -1, scal = FALSE)
}

\arguments{
  \item{x}{GVAR object}
  \item{n}{scalar, number of periods}
  \item{shock.var}{list of vectors of length 2, first element specifies the country to be shocked, second element the variable in the country (e.g. \code{shock.var=c(13,3)} third variable in US is going to be shocked)}
  \item{shock.dir}{a list of scalars, set to 1 or -1, for a (negative) 1sd shock, where sd refers to the country specific residuals. Alternatively, if scal = TRUE, any number submitted indicates the magnitude of the resulting first period shock, e.g. -0.01 for a 1\% decrease in the first quarter}
  \item{scal}{TRUE if shock.dir indicates resulting value of shock}
}
\details{ xxx
}
\value{ xxx
}
\references{xxx}
\author{
xxx
}
\note{
xxx
}


\seealso{\code{\link{GVAR}}}
\examples{
data(pesaran26)
c.names <- names(Data)[-length(Data)]

p <- c(2,2,2,1,2,2,1,2,2,2,2,1,2,1,1,2,2,2,2,2,2,1,2,2,2,2)
q <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
lex <- q

endo <- ord <- we <- d <- vector("list",length=length(c.names))
names(endo) <- names(ord) <- names(we) <- names(d) <- c.names
# base country: usa
endo[[1]] <- c(1:3,5:7)
ord[[1]] <- c(1:3,5:7)
we[[1]] <- c(1:2,4)
d[[1]] <- NULL
# countries with 6 endogenous variables:
for (j in c("EuroArea", "Japan", "UK", "Sweden", "Switzerland", "Norway", "Australia", "Canada", "NewZealand", "Korea", "Safrica")) 
{i <- which(c.names==j); endo[[i]] <- ord[[i]] <- 1:6}
# countries with 5 endogenous variables:
for (j in c("Argentina", "Chile", "Malaysia", "Philippines", "Singapore", "Thailand", "India")) 
{i <- which(c.names==j); endo[[i]] <- ord[[i]] <- 1:5}
# countries with 4 endogenous variables:
for (j in c("China", "Brazil", "Mexico", "Peru", "Indonesia", "Turkey")) 
{i <- which(c.names==j); endo[[i]] <- ord[[i]] <- c(1:2,4:5)}
# Saudi Arabia
endo[[21]] <- ord[[21]] <- c(1:2,4)

# all countries but us
for (i in 2:length(we))
{
  we[[i]] <- c(1:3,5,6)
  d[[i]] <- 1
}

Case <- "IV"
r <- c(2,1,1,4,3,3,3,2,2,1,2,3,3,4,4,3,3,4,1,2,3,3,2,1,1,1)

res.GVAR <- GVAR(Data=Data,r=r,p=p,q=q,weight=weight,Case=Case,exo.var=TRUE,d=d,lex=lex,ord=ord,we=we,endo=endo,method="max.eigen")

# OIR
testOIR <- OIR(res.GVAR,n=40,shock.var=1,shock.dir=-1,scal=FALSE)
}
\keyword{ts}

