est.we.mdls2 <- function (z.ts,etw,p,q=p,n,ex=0,lex=NULL,case,r,we.test=FALSE)
## notation: mixture of
# [1] Pesaran et al. "Structural Analysis of VECMs with exog. I(1) variables", Journal of Econometrics 97 (2000) p. 293-343
# [2] Johansen's book "Likelihood-based Inference in Cointegrated VAR Models"
## model for z_t = (y_t',x_t')':
# \Delta y_t = c_0+c_1*t+\Lambda\Delta x_t+\sum_{i=1}^{p-1}\Psi_i\Delta z_{t-i}+\Pi_y z_{t-1}+\u_t, where u_t is N(0,\Omega_{uu}), and
# \Delta x_t = a_{x0}+\sum_{i=1}^{p-1}\Gamma_{xi}\Delta z_{t-i}+e_{xt}, where e_{xt} is N(0,\Omega_{xx})
## to do/features:
# no include seasonal dummies
# ok estimate quarterly data
# ok include case p=1
# ok allow for higher AR-lag-order in \Delta x_t: up to \Delta x_{t-q+1}
{
## critical values are taken from [1]:
#source("cv.tables.pesaran.r")

freq<- etw[["freq"]] # time sampling frequency
dt<- 1/freq # time sampling interval
T<- (etw[["end"]]-etw[["start"]])*freq+1  # number of time samples for estimation
m<- dim(z.ts)[2]-ex # number of variables in z_t
# n is the number of endogenous variables y_t
k<- m-n # number of weakly exogenous variables x_t

if (is.null(lex)) lex <- 0

## rownames of time series y, x and z
if (is.null(dimnames(z.ts)[[2]])) {
 dimnames(z.ts)[[2]]<- c(paste("y",1:n,sep=""),paste("x",1:k,sep=""),paste("d",1:ex,sep=""))
}
y.ts<- z.ts[,1:n]
dimnames(y.ts)[[2]]<- dimnames(z.ts)[[2]][1:n]
x.ts <- z.ts[,(n+1):(n+k)]
if (k>1) {dimnames(x.ts)[[2]]<- dimnames(z.ts)[[2]][(n+1):(n+k)]}
d.ts <- z.ts[,-(1:m)]
if (ex>1) {dimnames(d.ts)[[2]]<- dimnames(z.ts)[[2]][-(1:m)]}

## construct deterministic terms (Dt) matrix: cf. [2]
if (case=="I") { # c.0=0 and c.1=0
 Dt<- rbind();
} else if (case=="II") { # c.0=-Pi.y\mu and c.1=0
} else if (case=="III") { # c.0!=0 and c.1=0
 Dt<- rbind(rep(1,T))
 rownames(Dt)<- c("constant")
} else if (case=="IV") { # c.0!=0 and c.1=-Pi.y\gamma
} else if (case=="V") { # c.0!=0 and c.1!=0
 Dt<- rbind(rep(1,T),1:T)
 rownames(Dt)<- c("constant","t")
}

# seq(etw[["start"]],etw[["end"]],by=dt)

## construct Z-matrices from data: cf. [2]
Z0<- t(diff(window(y.ts,start= etw[["start"]]-dt, end= etw[["end"]]))) # Delta y_t
rownames(Z0)<- paste("D",dimnames(y.ts)[[2]],"-0",sep="")
Z1<- t(window(z.ts[,1:m],start= etw[["start"]]-dt, end= etw[["end"]]-dt))    # z_{t-dt}
rownames(Z1)<- paste(rownames(Z1),"-1",sep="")
Z2<- t(diff(window(x.ts,start= etw[["start"]]-dt, end= etw[["end"]]))) # Delta x_t
rownames(Z2)[1:k]<- paste("D",colnames(z.ts)[(n+1):(n+k)],"-0",sep="") # paste("Dx",1:k,"-0",sep="")
if (max(p,q)>1) {
  if (min(p,q)>1) {
    for (i in 1:(min(p,q)-1)) { # include p-1 lags Delta z_{t-i}
      Z2 <- rbind(Z2, t(diff(window(z.ts[,1:m],start= etw[["start"]]-(1+i)*dt,end= etw[["end"]]-i*dt))))
      rownames(Z2)[(k+(i-1)*m+1):(k+i*m)] <- paste("D",dimnames(z.ts)[[2]][1:m],"-",i,sep="")  # paste("Dz", 1:m,"-",i, sep = "")
    }
  }
 if (p>q) {
  for (i in q:(p-1)) {
   Z2 <- rbind(Z2, t(diff(window(y.ts,start= etw[["start"]]-(1+i)*dt,end= etw[["end"]]-i*dt))))
   rownames(Z2)[(k+(min(p,q)-1)*m+(i-q)*n+1):(k+(min(p,q)-1)*m+(i-q)*n+n)]<- paste("D",dimnames(y.ts)[[2]],"-",i,sep="")
  }
 } else if (q>p) {
  for (i in p:(q-1)) {
   Z2 <- rbind(Z2, t(diff(window(x.ts,start= etw[["start"]]-(1+i)*dt,end= etw[["end"]]-i*dt))))
   rownames(Z2)[(k+(min(p,q)-1)*m+(i-p)*k+1):(k+(min(p,q)-1)*m+(i-p)*k+k)]<- paste("D",colnames(z.ts)[(n+1):(n+k)],"-",i,sep="")
  }
 }
}

if (ex!=0)
{
  temp <- t(window(d.ts,start= etw[["start"]]-dt, end= etw[["end"]]-dt))
  rownames(temp) <- paste(colnames(z.ts)[-(1:m)],"-1",sep="")
  Z1 <- rbind(Z1,temp)

  for (i in 0:(lex-1))
  {
    Z2 <- rbind(Z2,t(diff(window(d.ts,start=etw[["start"]]-(1+i)*dt,end=etw[["end"]]-i*dt))))
    rownames(Z2)[((k+(p-1)*n+(q-1)*k)+(i+1)*1):((k+(p-1)*n+(q-1)*k)+ex*(i+1))] <- paste("D",colnames(z.ts)[-(1:m)],"-",i,sep="")
  }
}


if (is.element(case, c("I","III","V"))) {
 Z2<- rbind(Z2, Dt)
} else if (case=="II") {
 Z1<- rbind(Z1, 1); rownames(Z1)[dim(Z1)[1]]<- "constant"
} else if (case=="IV") {
 Z1<- rbind(Z1,1:T)
 rownames(Z1)[dim(Z1)[1]]<- "t"
 Z2<- rbind(Z2, 1); rownames(Z2)[(p-1)*n+q*k+lex*ex+1]<- "constant"
} else {stop("\nUnkown case.\n")}

# if (!(is.null(season))) { # seasonal dummies: still to implement?
#  l<- (season.start.time-etw[["start"]])%%season
#  dum<- (diag(season) - 1/season)[-season,]
#  dum<- matrix(dum,nrow=nrow(dum),ncol=season*(ceiling(T/season)+1))
#  dum<- dum[,-(1:(season-l))]
#  dum<- dum[,1:T]
#  Z2<- rbind(Z2,dum)
# }

## product moment matrices M_{ij}, Residuals R_i: cf. [2]
M00 <- tcrossprod(Z0)/T # tcrossprod(x,y) is the same as x%*%t(y) but faster
M11 <- tcrossprod(Z1)/T
M22 <- tcrossprod(Z2)/T
M01 <- tcrossprod(Z0, Z1)/T
M02 <- tcrossprod(Z0, Z2)/T
M20 <- tcrossprod(Z2, Z0)/T
M10 <- tcrossprod(Z1, Z0)/T
M12 <- tcrossprod(Z1, Z2)/T
M21 <- tcrossprod(Z2, Z1)/T
# M11.inv <- solve(M11)
M22.inv <- solve(M22)
M22.inv <- (M22.inv+t(M22.inv))/2 # make M22.inv symmetric
R0 <- Z0-M02%*%M22.inv%*%Z2
R1 <- Z1-M12%*%M22.inv%*%Z2
S00 <- tcrossprod(R0,R0)/T
S01 <- tcrossprod(R0,R1)/T
S10 <- tcrossprod(R1,R0)/T
S11 <- tcrossprod(R1,R1)/T
S00.inv <- solve(S00)
S00.inv <- (S00.inv+t(S00.inv))/2 # make S00.inv symmetric
S11.inv <- solve(S11)
S11.inv <- (S11.inv+t(S11.inv))/2 # make S11.inv symmetric
## generalized eigenvalue problem: cf. [2]
# |\lambda S11-S10 S00.inv S01|=0
N<- S11
M<- S10%*%S00.inv%*%S01
M<- (M+t(M))/2 # make M symmetric

C <- t(chol(N)) # C%*%t(C)=N
C.inv<- solve(C)
eig<- eigen(C.inv%*%M%*%t(C.inv),symmetric=TRUE)
lambda<- eig[["values"]] # already sorted in decreasing order
V<- t(C.inv)%*%eig[["vectors"]]
## eventually normalize V:
V.orig <- V
if (0) {V <- sapply(1:m, function(x) V[,x]/V[1,x])}

if (r==0){
  beta <- NULL
  alpha <- NULL
  Pi.y <- matrix(0,nrow=n,ncol=dim(Z1)[1])
  Psi_ <- M02%*%M22.inv
  Omega.uu <- S00
}

if ((r>0) & (r<n)) {
 beta <- as.matrix(V[,1:r])
 print(beta)
 alpha <- S01%*%beta%*%solve(t(beta)%*%S11%*%beta)
 print(alpha)
 print(alpha%*%t(beta))
 if (1) {beta <- as.matrix(V[,1:r]%*%solve(V[1:r,1:r]))}
 print(beta)
 alpha <- S01%*%beta%*%solve(t(beta)%*%S11%*%beta)
 print(alpha)
 Pi.y <- alpha%*%t(beta)
 print(Pi.y)
 Psi_ <- M02%*%M22.inv-Pi.y%*%M12%*%M22.inv
 Omega.uu <-S00-Pi.y%*%S11%*%t(Pi.y)
}

if (r==n) {
  beta <- NULL
  alpha <- NULL
  Pi.y <- S01%*%S11.inv
  Psi_ <- M02%*%M22.inv-Pi.y%*%M12%*%M22.inv
  Omega.uu <-S00-Pi.y%*%S11%*%t(Pi.y)
}

## residuals ############################################################

U <- Z0 - Pi.y%*%Z1 - Psi_%*%Z2

# this is needed for t-values
Y <- Z0 - Pi.y%*%Z1 

## compute Psi_i, c.x, etc matrices and vectors
Lambda <- as.matrix(Psi_[,1:k])
Phi <- NULL
Psi <- NULL
if (max(q,lex)>1) {Psis <- vector("list", max(q,lex)-1)}
if (max(p,q)>1) {
 Phis <- vector("list", p-1)
 if (min(p,q)>1) {
  for (i in (1:(min(p,q)-1)) ) {
   Phis[[i]] <- as.matrix(Psi_[,(k+1+(i-1)*m):(k+(i-1)*m+n)])
   Psis[[i]] <- as.matrix(Psi_[,(k+1+n+(i-1)*m):(k+n+(i-1)*m+k)])
   colnames(Psis[[i]]) <- colnames(Psi_)[(k+1+n+(i-1)*m):(k+n+(i-1)*m+k)]
  }
 }
 if (p>q) {
  for (i in (q:(p-1)) ) {
   Phis[[i]] <- as.matrix(Psi_[,(k+1+(q-1)*m+(i-q)*n):(k+(q-1)*m+(i-q)*n+n)])
  }
 }
 if (q>p) {
  for (i in (p:(q-1)) ) {
   Psis[[i]] <- as.matrix(Psi_[,(k+1+(p-1)*m+(i-p)*k):(k+(p-1)*m+(i-p)*k+k)])
   colnames(Psis[[i]]) <- colnames(Psi_)[(k+1+(p-1)*m+(i-p)*k):(k+(p-1)*m+(i-p)*k+k)]
  }
 }
 if (p>1) {Phi <- Phis}
 if (q>1) {Psi <- Psis}
}

if (ex!=0)
{
  Lambda <- cbind(Lambda,as.matrix(Psi_[,((p-1)*n+q*k+1):((p-1)*n+q*k+ex)],ncol=ex))
  colnames(Lambda) <- colnames(Psi_)[c(1:k,((p-1)*n+q*k+1):((p-1)*n+q*k+ex))]
  if (lex>1)
  {
  	for (i in 1:(lex-1))
    {
      if (i<=(q-1))
      {
        Psi[[i]] <- cbind(Psi[[i]],as.matrix(Psi_[,((p-1)*n+q*k+i*ex+1):((p-1)*n+q*k+(i+1)*ex)],ncol=ex))
        colnames(Psi[[i]])[-(1:k)] <- colnames(Psi_)[((p-1)*n+q*k+i*ex+1):((p-1)*n+q*k+(i+1)*ex)]
      } else {
      	Psi[[i]] <- as.matrix(Psi_[,((p-1)*n+q*k+i*ex+1):((p-1)*n+q*k+(i+1)*ex)],ncol=ex)
      	colnames(Psi[[i]]) <- colnames(Psi_)[((p-1)*n+q*k+i*ex+1):((p-1)*n+q*k+(i+1)*ex)]
      }
    }
  }
}


c.0 <- NULL
c.1 <- NULL
if (case=="V") {
  c.0 <- as.matrix(Psi_[,dim(Psi_)[2]-1],n,1)
  colnames(c.0) <- "Const"
  c.1 <- as.matrix(Psi_[,dim(Psi_)[2]],n,1)
  colnames(c.1) <- "Trend"
} else if (case=="III") {
  c.0 <- as.matrix(Psi_[,dim(Psi_)[2]],n,1)
  colnames(c.0) <- "Const"
} else if (case=="I") {
} else if (case=="IV") {
  c.0 <- as.matrix(Psi_[,dim(Psi_)[2]],n,1)
  colnames(c.0) <- "Const"
  if (r>0 & r<n) {
    c.1 <- alpha%*%t(beta)[,m+ex+1]
    colnames(c.1) <- "Trend"
#    beta <- beta[-(m+1),]
  }
  c.1 <- as.matrix(Pi.y[,m+ex+1],n,1)
  colnames(c.1) <- "Trend"
#  Pi.y <- Pi.y[,1:m]
} else if (case=="II") {
  if (r>0 & r<n) {
    c.0 <- alpha%*%t(beta)[,m+ex+1]
    colnames(c.0) <- "Const"
#    beta <- beta[-(m+1),]
  }
  c.0 <- as.matrix(Pi.y[,m+ex+1],n,1)
  colnames(c.0) <- "Const"
#  Pi.y <- Pi.y[,1:m]
}

## test for weak exogeneity #############################################################

we.test.res <- NULL
if (we.test==TRUE)
{
  we.res <- matrix(NA,nrow=m-n,ncol=3)
  colnames(we.res) <- c("F Stat.","crit. Value","p-Value")
  rownames(we.res) <- colnames(z.ts)[(n+1):m]
  for (i in 1:(m-n))
  {
  	t1 <- lm(t(Z2)[,i] ~ t(Z2)[,-(1:(m-n))])
  	t2 <- lm(t(Z2)[,i] ~ t(t(beta)%*%Z1)+t(Z2)[,-(1:(m-n))])
  	
  	Fstat <- ((sum(t1$residuals^2)-sum(t2$residuals^2))/dim(t(t(beta)%*%Z1))[2])/(sum(t2$residuals^2)/(T-dim(t(t(beta)%*%Z1))[2]-dim(t(Z2)[,-(1:(m-n))])[2]))
  	critV <- qf(0.95,dim(t(t(beta)%*%Z1))[2],T-dim(t(t(beta)%*%Z1))[2]-dim(t(Z2)[,-(1:(m-n))])[2])
  	pValue <- 1-pf(Fstat,dim(t(t(beta)%*%Z1))[2],T-dim(t(t(beta)%*%Z1))[2]-dim(t(Z2)[,-(1:(m-n))])[2])
  	cat("F(",dim(t(t(beta)%*%Z1))[2],",",T-dim(t(t(beta)%*%Z1))[2]-dim(t(Z2)[,-(1:(m-n))])[2],")\n")
	
	 we.res[i,] <- c(Fstat,critV,pValue)
  }
  we.test.res <- we.res
}

## t-values #############################################################

se <- new.env()
tvals <- new.env()
pvals <- new.env()

# beta
if (!(r==n))
{
  M <- diag(T)-t(Z2)%*%solve(Z2%*%t(Z2))%*%Z2
  Sigma.u.tilde <- U%*%t(U)/T
  Omega.b <- solve(Z1[-(1:r),]%*%M%*%t(Z1[-(1:r),]))%x%solve(t(alpha)%*%solve(Sigma.u.tilde)%*%alpha)
  tvals$beta <- beta[-(1:r),]/sqrt(diag(Omega.b))
}

# other parameters

res <- lm(t(Y) ~ t(Z2) + 0)
weights <- log((res$residuals)^2)

res.w <-  lm(weights ~ t(Z2) + 0)
weights.final <- 1/(exp(res.w$fitted.values))
weights <- data.frame(weights.final)

res.final <- lm(t(Y) ~ t(Z2) + 0, weights)

if (k>=1)
{
  if (ex!=0)
  {
    se$Lambda <- matrix(NA,n,k+ex)
    tvals$Lambda <- matrix(NA,n,k+ex)
    pvals$Lambda <- matrix(NA,n,k+ex)
    for (i in 1:n)
    {
      se$Lambda[i,] <- summary(res.final)[[i]]$coefficients[c(1:k,((p-1)*n+q*k+1):((p-1)*n+q*k+ex)),2]
      tvals$Lambda[i,] <- summary(res.final)[[i]]$coefficients[c(1:k,((p-1)*n+q*k+1):((p-1)*n+q*k+ex)),3]
      pvals$Lambda[i,] <- summary(res.final)[[i]]$coefficients[c(1:k,((p-1)*n+q*k+1):((p-1)*n+q*k+ex)),4]      
    }  
  } else {
    se$Lambda <- matrix(NA,n,k)
    tvals$Lambda <- matrix(NA,n,k)
    pvals$Lambda <- matrix(NA,n,k)
    for (i in 1:n)
    {
      se$Lambda[i,] <- summary(res.final)[[i]]$coefficients[1:k,2]
      tvals$Lambda[i,] <- summary(res.final)[[i]]$coefficients[1:k,3]
      pvals$Lambda[i,] <- summary(res.final)[[i]]$coefficients[1:k,4]      
    }
  }
  dimnames(se$Lambda) <- dimnames(tvals$Lambda) <- dimnames(pvals$Lambda) <- dimnames(Lambda)
}

if (p>1)
{
  se$Phi <- list()
  tvals$Phi <- list()
  pvals$Phi <- list()
  qp <- 1
  for (j in 1:(p-1)) 
  {
    se$Phi[[j]] <- matrix(NA,n,n)
    tvals$Phi[[j]] <- matrix(NA,n,n)
    pvals$Phi[[j]] <- matrix(NA,n,n)
    for (i in 1:n)
    {
      
      se$Phi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+qp+1):(j*n+qp),2]
      tvals$Phi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+qp+1):(j*n+qp),3]
      pvals$Phi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+qp+1):(j*n+qp),4]      
    }
    dimnames(se$Phi[[j]]) <- dimnames(tvals$Phi[[j]]) <- dimnames(pvals$Phi[[j]]) <- dimnames(Phi[[j]])
    if ((q-1)>=j) {qp <- qp+1}
  }
}

if (max(q,lex)>1)
{
  se$Psi <- list()
  tvals$Psi <- list()
  pvals$Psi <- list()
  qp <- 1
  for (j in 1:(max(q,lex)-1)) 
  {
    se$Psi[[j]] <- matrix(NA,n,dim(Psi[[j]])[2])
    tvals$Psi[[j]] <- matrix(NA,n,dim(Psi[[j]])[2])
    pvals$Psi[[j]] <- matrix(NA,n,dim(Psi[[j]])[2])
    for (i in 1:n)
    {
      if ((min(q,lex)-1>=1) && ex!=0)
      {
        index <- c((qp*n+1+k*j):(qp*n+k*(j+1)),((p-1)*n+q*k+ex*j+1):((p-1)*n+q*k+ex*j+ex))
      } else if (q>lex) {
      	index <- (qp*n+1+k*j):(qp*n+k*(j+1))
      } else {
        index <- ((p-1)*n+q*k+ex*j+1):((p-1)*n+q*k+ex*j+ex)
      }
      se$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[index,2]
      tvals$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[index,3]
      pvals$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[index,4]      
    }
    dimnames(se$Psi[[j]]) <- dimnames(tvals$Psi[[j]]) <- dimnames(pvals$Psi[[j]]) <- dimnames(Psi[[j]])
    if ((p-1)>=j) {qp <- qp+1}
  }
}

if (case=="III" || case=="IV" || case=="V")
{
  temp <- 1 
  se$c.0 <- matrix(NA,n,1)    
  tvals$c.0 <- matrix(NA,n,1)
  pvals$c.0 <- matrix(NA,n,1)                         
  for (i in 1:n)
  {
    se$c.0[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+k*q+lex*ex+temp),2]
    tvals$c.0[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+k*q+lex*ex+temp),3]
    pvals$c.0[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+k*q+lex*ex+temp),4]
  }
  temp <- temp+1
} 
if (case=="V") 
{
  se$c.1 <- matrix(NA,n,1)    
  tvals$c.1 <- matrix(NA,n,1)
  pvals$c.1 <- matrix(NA,n,1)
  for (i in 1:n)
  {
    se$c.1[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+(m-n)*q+lex*ex+temp),2]
    tvals$c.1[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+(m-n)*q+lex*ex+temp),3]
    pvals$c.1[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+(m-n)*q+lex*ex+temp),4]
  }
  temp <- temp+1
}

#if (!is.null(season))
#{
#  se$season <- matrix(NA,n,season-1)    
#  tvals$season <- matrix(NA,n,season-1)
#  pvals$season <- matrix(NA,n,season-1)
#  for (i in 1:n)
#  { 
#  se$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+temp):dim(summary(res.final)[[i]]$coefficients)[1],2]
#  tvals$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+temp):dim(summary(res.final)[[i]]$coefficients)[1],3]
#  pvals$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+temp):dim(summary(res.final)[[i]]$coefficients)[1],4]
#  }
#  rownames(se$season) <- rownames(tvals$season) <- rownames(pvals$season) <- rownames(Phi)
#  colnames(se$season) <- colnames(tvals$season) <- colnames(pvals$season) <- paste("s",1:(season-1),sep="")
#}
se <- as.list(se)
tvals <- as.list(tvals)
pvals <- as.list(pvals)

## output:
mdls <- list(type="weakly exogenous VECM",dat=z.ts,freq=freq,m=m,n=n,p=p,q=q,ex=ex,lex=lex,r=r,T=T,alpha=alpha,beta=beta,Pi.y=Pi.y,Phi=Phi,Psi=Psi,Lambda=Lambda,case=case,c.0=c.0,c.1=c.1,Omega.uu=Omega.uu,S=list(S00=S00,S10=S10,S01=S01,S11=S11),lambda=lambda,residuals=t(U),se=se,tvals=tvals,pvals=pvals,we.test.res=we.test.res)
class(mdls) <- "vecm"
return(mdls)

}

