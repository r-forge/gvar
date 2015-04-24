\name{GVAR}
\alias{GVAR}
\title{Global Vector Auto-Regressive Modelling}
\description{
  \code{GVAR} computes VECMs for all regions and stacks the models to a Global Vector Autoregressive Model
}
\usage{
GVAR(Data, tw = NULL, p, q = p, r = NULL, weight, Case, exo.var = FALSE,
     d = NULL, lex = NULL, endo = NULL, ord = NULL, we = NULL, method = "max.eigen", caseTest = FALSE, weTest = FALSE)
}

\arguments{
  \item{Data}{timeseries data as list (each entry is a matrix of a subsystem of variables, if \code{exo.var=TRUE} the last entry are exogeneous variables)}
  \item{tw}{time window, vector of start and end point, if \code{NULL} the maximum time interval will be used}
  \item{p}{scalar/vector of endogenous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{q}{scalar/vector of (weakly) exogeneous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{r}{scalar/vector of cointegrating relations, if a scalar is provided the same cointegration rank is used for all subsystems, if set to \code{NULL} the cointegration rank will be determind by \code{method}}
  \item{weight}{weight matrix, the diagonal elements need to be zero}
  \item{Case}{single value/vector of cases (\code{I} to \code{V)}, where case \code{I} is a zero intercept, zero trend model, case \code{II} is a restricted intercept, zero trend model, \code{III} is a unrestricted intercept, zero trend model, \code{IV} is a unrestricted intercept restricted trend model and \code{V} is a unrestricted intercept, unrestricted trend model; if a single value is provided the same structure is used for all subsystems}
  \item{exo.var}{if \code{TRUE} strictly exogeneous variables are included in the model}
  \item{d}{list showing which strictly exogeneous variables enter the subsystem equations, if \code{NULL} all variables in the last entry of \code{data} will be used for all subsystems}
  \item{lex}{scalar/vector of strictly exogeneous lags, if a scalar is provided the same lag length is used for all subsystems}
  \item{endo}{list of endogenous variables used in each subsystem, if \code{NULL} all variables in \code{data} will be used}
  \item{ord}{vector used if variables in the different subsystem don't appear in the same order, order of each subsystem is concatenated to one vector, if \code{NULL} the variables in \code{data} are assumed to be ordered identically in all subsystems}
  \item{we}{list with numbers of weakly exogeneous variables included in each VECM, corresponds to numbers in \code{ord}, if \code{NULL} every variable appearing in all subsystems will be used}
  \item{method}{select cointegrating rank by max. eigenvalue (\code{max.eigen}) or trace statistic (\code{trace})}
  \item{caseTest}{provide test statistics regarding the intercept/trend structure}
  \item{weTest}{perform F test for weak exogenity}
}
\details{
  The function computes a VECM for every subsystem before stacking the results to a GVAR model.
  
  Specification of input here.
}
\value{
  An object of class \code{GVAR} containing the following items:
  
  \item{subsys}{subsystem names}
  \item{Data}{data}
  \item{we.vecms}{VECMs of the subsystems}
  \item{X}{data as one single matrix}
  \item{bigT}{length of time series data}
  \item{r}{vector of cointegration ranks of the VECMs}
  \item{Case}{vector of intercept/trend behaviour of the VECMs}
  \item{W}{multiplier matrix to generate endogenous and weakly exogenous variables from X}
  \item{G}{multiplier matrix for the current variables}
  \item{H}{multiplier matrix for the lagged variables}
  \item{Upsilon.0}{multiplier matrix for the current strictly exogenous variables}
  \item{Upsilon}{multiplier matrix for the lagged strictly exogenous variables}
  \item{c.0}{multiplier matrix for the intercept}
  \item{c.1}{multiplier matrix for the trend}
  \item{caseTest}{test statistics for case selection}
  \item{weight}{weight matrix used to calculate the weakly exogenous variables}
  \item{U}{residuals of the GVAR}
  \item{U.cov}{residual covariance matrix}
  \item{arguments}{arguments passed to GVAR function, including lags, variable types,...}
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

\seealso{\code{\link{est.we.mdls}}}

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

# view vecm models
res.GVAR$we.vecms
}

\keyword{ts}

