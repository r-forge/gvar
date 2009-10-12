## fctns.spec.tests.r

ur.tests <- function(X)
{

tests<- list()

## use package tseries:
library(tseries)
if (!is.matrix(X)){X <- as.matrix(X,ncol=1)}

adf <- list()
pp <- list()
kpss <- list()

for (i in 1:dim(X)[2])
{
  adf[[i]] <- adf.test(X[,i])  # Augmented Dickey-Fuller Test with null hypothesis that x has a unit root
  pp[[i]] <-pp.test(X[,i])     # Phillips-Perron Unit Root Test with null hypothesis that x has a unit root
  kpss[[i]] <-kpss.test(X[,i]) # Kwiatkowski-Phillips-Schmidt-Shin (KPSS) Test with null hypothesis that x is level or trend stationary
}

kpss[[1]]$alternative <- "unit root"

tests$all<-list(adf=adf,pp=pp,kpss=kpss)

essence <- array(dim=c(length(tests$all),dim(X)[2]+1))

values <- list()
nam <- vector()
for (k in 1:length(tests$all))
{
  values[[k]] <- vector()
  for (l in 1:dim(X)[2])
  {
    values[[k]][l] <- round(tests$all[[k]][[l]]$p.value,digits=5)
  }
  essence[k,]<- c(tests$all[[k]][[1]]$alternative, values[[k]])
  nam[k] <- tests$all[[k]][[1]]$method
}

essence<- as.data.frame(essence)
names(essence)<- c("  alternative",paste("  p-value",paste("Var.",1:dim(X)[2],sep="")))
dimnames(essence)[[1]] <- nam

return(essence)
cat("\n")
}