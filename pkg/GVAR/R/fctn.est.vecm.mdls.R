est.vecm.mdls <- function (Y.ts,etw,p,case,r,ex=0,lex=NULL,season=NULL,season.start.time=NULL)
## notation: almost as in Johansen's book "Likelihood-based Inference in Cointegrated VAR Models"
## model:
# \Delta Y_t = \Pi Y_{t-1}+\sum_{i=1}^{k-1} \Gamma_i\Delta Y_{t-i}+\mu_0+\mu_1 t+\Phi D_t+\epsilon_t, where \epsilon_t is N(0,\Omega)
{
freq <- etw[["freq"]] # time sampling frequency
dt <- 1/freq # time sampling interval
T <- (etw[["end"]]-etw[["start"]])*freq+1  # number of time samples for estimation
n <- dim(Y.ts)[2] - ex # number of endogenous variables in Y_t

if (is.null(lex)) lex <- 0

if (ex>0)
{
  d.ts <- Y.ts[,-(1:n)]
  Y.ts <- Y.ts[,1:n]
  if (ex>1) {dimnames(d.ts)[[2]]<- dimnames(Y.ts)[[2]][-(1:n)]}
}

if (case=="I"){case <- "H_2(r)"}
if (case=="II"){case <- "H_1^*(r)"}
if (case=="III"){case <- "H_1(r)"}
if (case=="IV"){case <- "H^*(r)"}
if (case=="V"){case <- "H(r)"}

## construct deterministic terms (Dt) matrix:
if (case=="H(r)") {
 Dt<- rbind(rep(1,T),seq(etw[["start"]],etw[["end"]],by=dt))
 rownames(Dt)<- c("Const","t")                             
} else if (case=="H_1(r)") {
 Dt<- rbind(rep(1,T)); rownames(Dt)<- c("Const")
} else if (case=="H_2(r)") {
 Dt<- rbind();
} else if (case=="H^*(r)") {
 Dt<- rbind(rep(1,T),seq(etw[["start"]],etw[["end"]],by=dt))
 rownames(Dt)<- c("Const","t")
} else if (case=="H_1^*(r)") {
 Dt<- rbind(rep(1,T)); rownames(Dt)<- c("Const")
}
## construct Z-matrices from data
# Y<- t(window(Y.ts,start = etw[["start"]], end = etw[["end"]]))
Z0<- t(diff(window(Y.ts,start= etw[["start"]]-dt, end= etw[["end"]])))
Z1<- t(window(Y.ts,start= etw[["start"]]-dt, end= etw[["end"]]-dt))
Z2<- matrix(NA, nrow=0, ncol=T)

if (p>1) {
 for (i in 1:(p-1)) {
  Z2<- rbind(Z2, t(diff(window(Y.ts,start= etw[["start"]]-(1+i)*dt , end= etw[["end"]]-i*dt))))
  rownames(Z2)[((i-1)*n+1):(i*n)]<- paste("DY", 1:n,"-",i, sep = "")
 }
}

if (ex!=0)
{
  for (i in 0:(lex-1))
  {
    Z2 <- rbind(Z2,t(diff(window(d.ts,start=etw[["start"]]-(1+i)*dt,end=etw[["end"]]-i*dt))))
    rownames(Z2)[((p-1)*n+ex*i+1):((p-1)*n+ex*(i+1))] <- paste("D",colnames(Y.ts)[-(1:n)],"-",i,sep="")
  }
}


if (is.element(case, c("H(r)","H_1(r)","H_2(r)"))) {
 Z2<- rbind(Z2, Dt)
} else if (case=="H^*(r)") {
 Z1<- rbind(Z1, seq(etw[["start"]],etw[["end"]],by=dt))
 rownames(Z1)[n+1]<- "t"
 Z2<- rbind(Z2, 1)
 rownames(Z2)[(p-1)*n+1]<- "Const"
} else if (case=="H_1^*(r)") {
 Z1<- rbind(Z1, 1)
 rownames(Z1)[n+1]<- "Const"
} else {stop("\nUnkown case.\n")}

if (!(is.null(season))) { # seasonal dummies
 l<- (season.start.time*freq-etw[["start"]]*freq)%%season
 dum<- diag(season)[-season,]
 # dum<- (diag(season) - 1/season)[-season,]
 dum<- matrix(dum,nrow=nrow(dum),ncol=season*(ceiling(T/season)+1))
 dum<- dum[,-(1:(season-l))]
 dum<- dum[,1:T]
 Z2<- rbind(Z2,dum)
}
## product moment matrices M_{ij}, Residuals R_i
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
if (length(Z2)) {
 M22.inv <- solve(M22)
 R0 <- Z0-M02%*%M22.inv%*%Z2
 R1 <- Z1-M12%*%M22.inv%*%Z2
} else {
 R0 <- Z0
 R1 <- Z1
}
S00 <- tcrossprod(R0,R0)/T
S01 <- tcrossprod(R0,R1)/T
S10 <- tcrossprod(R1,R0)/T
S11 <- tcrossprod(R1,R1)/T
S00.inv <- solve(S00)
S11.inv <- solve(S11)
## generalized eigenvalue problem
# |\lambda S11-S10 S00.inv S01|=0
N<- S11
M<- S10%*%S00.inv%*%S01
C <- t(chol(N)) # C%*%t(C)=N
C.inv<- solve(C)
eig<- eigen(C.inv%*%M%*%t(C.inv))
lambda<- eig[["values"]] # already sorted in decreasing order
V<- t(C.inv)%*%eig[["vectors"]]
## eventually normalize V
V.orig <- V
if (0) {V <- sapply(1:(n+ex), function(x) V[,x]/V[1,x])}

## compute beta, alpha, Pi, Psi, Omega etc. matrices for all ranks 0:p
if (r==0) {
 beta <- NULL
 alpha <- NULL
 Pi <- matrix(0,nrow=n,ncol=dim(Z1)[1])
 if (length(Z2)) Psi_ <- M02%*%M22.inv
 Omega <- S00
}

if ((r>0)&(r<n+ex)) {
 beta <- as.matrix(V[,1:r])
 if (1) {
 	beta <- as.matrix(V[,1:r]%*%solve(V[1:r,1:r]),ncol=r)
 	rownames(beta) <- rownames(V)
 	}
 alpha <- S01%*%beta%*%solve(t(beta)%*%S11%*%beta)
 Pi <- alpha%*%t(beta)
 if (length(Z2)) Psi_ <- M02%*%M22.inv-Pi%*%M12%*%M22.inv
 Omega <- S00-Pi%*%S11%*%t(Pi)
}

if (r==n+ex) {                                             
 beta <- NULL
 alpha <- NULL
 Pi <- S01%*%S11.inv
 if (length(Z2)) Psi_ <- M02%*%M22.inv-Pi%*%M12%*%M22.inv
 Omega <- S00-Pi%*%S11%*%t(Pi)
}

## residuals ############################################################

if (length(Z2)) {U <- Z0 - Pi%*%Z1 - Psi_%*%Z2} else {U <- Z0 - Pi%*%Z1}

# this is needed for t-values
Y <- Z0 - Pi%*%Z1 

## compute Gamma_i, mu_t, Psi, Phi matrices and vectors 

Gamma <- NULL
if (p>1) {
 Gamma_ <- list()
 for (i in (1:(p-1)) ) {
  Gamma_[[i]]<- Psi_[,(1+(i-1)*n):(i*n)]
 }
 Gamma <- Gamma_
}

mu0 <- NULL       # constant
mu1 <- NULL       # trend
Phi <- NULL       # season
Psi <- NULL       # exogenous

if ( case=="H(r)" && length(Z2) ) {
 if (ex>0) Psi <- Psi_[,((p-1)*n+1):((p-1)*n+ex*lex)]
 mu0 <- Psi_[,(p-1)*n+(lex-1)*ex+1]
 mu1 <- Psi_[,(p-1)*n+(lex-1)*ex+2]
 Phi <- Psi_[,-(1:((p-1)*n+lex*ex+2))]
 
} else if ( case=="H_1(r)" && length(Z2) ) {
 if (ex >0) Psi <- Psi_[,((p-1)*n+1):((p-1)*n+ex*lex)]
 mu0 <- Psi_[,(p-1)*n+(lex-1)*ex+1]
 Phi <- Psi_[,-(1:((p-1)*n+lex*ex+1))]
 
} else if ( case=="H_2(r)" && length(Z2) ) {
	if (p>1 || lex>0) 
  {
    if (ex>0) Psi <- Psi_[,((p-1)*n+1):((p-1)*n+ex*lex)]
    Phi <- Psi_[,-(1:((p-1)*n)+lex*ex)]
  } else {    
    Phi <- Psi_[,-(1:((p-1)*n)+lex*ex)]
  }
  
} else if ( case=="H^*(r)" ) {
 if ( length(Z2) ) {
  if (ex >0) Psi <- Psi_[,((p-1)*n+1):((p-1)*n+ex*lex)]
  mu0 <- Psi_[,(p-1)*n+lex*ex+1]
  Phi <- Psi_[,-(1:((p-1)*n+lex*ex+1))]
 }
 if (r>0 & r<n) {
  mu1 <- alpha%*%t(beta)[,n+1]
 }

} else if (case=="H_1^*(r)") {
 if ( length(Z2) ) {
  if (p>1 || lex>0) 
  {
    if (ex>0) Psi <- Psi_[,((p-1)*n+1):((p-1)*n+ex*lex)]
    Phi <- Psi_[,-(1:((p-1)*n+lex*ex))]
  } else {
    Phi <- Psi_
  }
 }
 if (r>0 & r<n) {
  mu0 <- alpha%*%t(beta)[,n+1]
 }
}

# build list of Psi

if (!is.null(Psi)) 
{
  Psi <- matrix(Psi,nrow=n)
  Phi_ <- vector("list",length=lex)
  for (i in 1:lex)
  {
    Phi_[[i]] <- matrix(Psi[,((i-1)*ex+1):(i*ex)],nrow=n)
    colnames(Phi_[[i]]) <- colnames(Psi_)[((p-1)*n+1+(i-1)*ex):((p-1)*n+ex*i)]
    rownames(Phi_[[i]]) <- rownames(Psi_) 
  }
  Psi <- Phi_
}

## t-values #############################################################

se <- new.env()
tvals <- new.env()
pvals <- new.env()

# beta
if (length(Z2)) {M <- diag(T)-t(Z2)%*%solve(Z2%*%t(Z2))%*%Z2} else {M <- diag(T)}

Sigma.u.tilde <- U%*%t(U)/T
if (r==1)
{
  Omega.b <- solve(matrix(Z1[-1,],ncol=T)%*%M%*%t(matrix(Z1[-1,],ncol=T)))%x%solve(t(alpha)%*%solve(Sigma.u.tilde)%*%alpha)
} else {
  Omega.b <- solve(matrix(Z1[-(1:r),],ncol=T)%*%M%*%t(matrix(Z1[-(1:r),],ncol=T)))%x%solve(t(alpha)%*%solve(Sigma.u.tilde)%*%alpha)
}
tvals$beta <- beta[-(1:r),]/sqrt(diag(Omega.b))
se$beta <- matrix(sqrt(diag(Omega.b)),ncol=r)

# other parameters
if (length(Z2))
{
  res <- lm(t(Y) ~ t(Z2) + 0)
  weights <- log((res$residuals)^2)

  res.w <-  lm(weights ~ t(Z2) + 0)  
  weights.final <- 1/(exp(res.w$fitted.values))
  weights <- data.frame(weights.final)

  res.final <- lm(t(Y) ~ t(Z2) + 0, weights)
}

if (p>1)
{
  se$Gamma <- list()
  tvals$Gamma <- list()
  pvals$Gamma <- list()
  for (j in 1:(p-1)) 
  {
    se$Gamma[[j]] <- matrix(NA,n,n)
    tvals$Gamma[[j]] <- matrix(NA,n,n)
    pvals$Gamma[[j]] <- matrix(NA,n,n)
    for (i in 1:n)
    {
      se$Gamma[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+1):(j*n),2]
      tvals$Gamma[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+1):(j*n),3]
      pvals$Gamma[[j]][i,] <- summary(res.final)[[i]]$coefficients[((j-1)*n+1):(j*n),4]      
    }
    dimnames(se$Gamma[[j]]) <- dimnames(tvals$Gamma[[j]]) <- dimnames(pvals$Gamma[[j]]) <- dimnames(Gamma[[j]])
  }
}

if (ex!=0)
{
  se$Psi <- list()
  tvals$Psi <- list()
  pvals$Psi <- list()
  for (j in 1:lex) 
  {
    se$Psi[[j]] <- matrix(NA,n,ex)
    tvals$Psi[[j]] <- matrix(NA,n,ex)
    pvals$Psi[[j]] <- matrix(NA,n,ex)
    for (i in 1:n)
    {
      se$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((p-1)*n+(j-1)*ex+1):((p-1)*n+j*ex),2]
      tvals$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((p-1)*n+(j-1)*ex+1):((p-1)*n+j*ex),3]
      pvals$Psi[[j]][i,] <- summary(res.final)[[i]]$coefficients[((p-1)*n+(j-1)*ex+1):((p-1)*n+j*ex),4]      
    }
    dimnames(se$Psi[[j]]) <- dimnames(tvals$Psi[[j]]) <- dimnames(pvals$Psi[[j]]) <- dimnames(Psi[[j]])
  }
}

temp <- 1
if (case=="H_1(r)" || case=="H^*(r)" || case=="H(r)")
{
  se$mu0 <- vector()    
  tvals$mu0 <- vector()
  pvals$mu0 <- vector()                         
  for (i in 1:n)
  {
    se$mu0[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),2]
    tvals$mu0[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),3]
    pvals$mu0[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),4]
  }
  temp <- temp+1
} 
if (case=="H(r)") 
{
  se$mu1 <- vector()    
  tvals$mu1 <- vector()
  pvals$mu1 <- vector()
  for (i in 1:n)
  {
    se$mu1[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),2]
    tvals$mu1[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),3]
    pvals$mu1[i] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp),4]
  }
  temp <- temp+1
}
if (!is.null(season))
{
  se$season <- matrix(NA,n,season-1)    
  tvals$season <- matrix(NA,n,season-1)
  pvals$season <- matrix(NA,n,season-1)
  for (i in 1:n)
  { 
  se$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp):dim(summary(res.final)[[i]]$coefficients)[1],2]
  tvals$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp):dim(summary(res.final)[[i]]$coefficients)[1],3]
  pvals$season[i,] <- summary(res.final)[[i]]$coefficients[(n*(p-1)+ex*lex+temp):dim(summary(res.final)[[i]]$coefficients)[1],4]
  }
  rownames(se$season) <- rownames(tvals$season) <- rownames(pvals$season) <- rownames(Phi)
  colnames(se$season) <- colnames(tvals$season) <- colnames(pvals$season) <- paste("s",1:(season-1),sep="")
}



se <- as.list(se)
tvals <- as.list(tvals)
pvals <- as.list(pvals)
 
                                                                                        
## output ###############################################################
# estimates of: alpha, beta, Pi, Gamma, mu_t, Omega
# t-values of beta, cf. p.100 in Luetkepohl
# various test statistics and critical values

# collect all estimated models for r=0:n in a list:
mdls <- list(type="pure VECM",dat=Y.ts,freq=freq,n=n,p=p,ex=ex,lex=lex,T=T,r=r,season=season,season.start.time=season.start.time,alpha=alpha,beta=beta,Pi=Pi,Gamma=Gamma,case=case,mu0=mu0,mu1=mu1,Phi=Phi,Psi=Psi,Omega=Omega,residuals=t(U),S=list(S00=S00,S10=S10,S01=S01,S11=S11),lambda=lambda,se=se,tvals=tvals,pvals=pvals) 
class(mdls) <- "vecm"
return(mdls)
}
