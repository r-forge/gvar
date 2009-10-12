rank.test.we <- function (z.ts,etw,p,q=p,n,ex,lex,case)
## notation: mixture of
# [1] Pesaran et al. "Structural Analysis of VECMs with exog. I(1) variables", Journal of Econometrics 97 (2000) p. 293-343
# [2] Johansen's book "Likelihood-based Inference in Cointegrated VAR Models"
## model for z_t = (y_t',x_t')':
# \Delta y_t = c_0+c_1*t+\Lambda\Delta x_t+\sum_{i=1}^{p-1}\Psi_i\Delta z_{t-i}+\Pi_y z_{t-1}+\u_t, where u_t is N(0,\Omega_{uu}), and
# \Delta x_t = a_{x0}+\sum_{i=1}^{p-1}\Gamma_{xi}\Delta z_{t-i}+e_{xt}, where e_{xt} is N(0,\Omega_{xx})
{
## critical values are taken from [1]:
data(cv.tables.pesaran)

freq<- etw[["freq"]] # time sampling frequency
dt<- 1/freq # time sampling interval
T<- (etw[["end"]]-etw[["start"]])*freq+1  # number of time samples for estimation
m<- dim(z.ts)[2]-ex # number of variables in z_t
# n is the number of endogenous variables y_t
k<- m-n # number of weakly exogenous variables x_t

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
   rownames(Z2)[(k+(min(p,q)-1)*m+(i-p)*k+1):(k+(min(p,q)-1)*m+(i-p)*k+k)]<- paste("D",colnames(z.ts)[-(c(1:n,(m+1):dim(z.ts)[2]))],"-",i,sep="")
  }
 }
}

if (ex!=0)
{
  for (i in 0:(lex-1))
  {
    Z2 <- rbind(Z2,t(diff(window(d.ts,start=etw[["start"]]-(1+i)*dt,end=etw[["end"]]-i*dt))))
    rownames(Z2)[((k+(p-1)*n+(q-1)*k)+(i+1)*1):((k+(p-1)*n+(q-1)*k)+ex*(i+1))] <- paste("D",colnames(z.ts)[-(1:m)],"-",i,sep="")
  }
}


if (is.element(case, c("I","III","V"))) {
 Z2<- rbind(Z2, Dt)
} else if (case=="II") {
 Z1<- rbind(Z1, 1); rownames(Z1)[m+1]<- "constant"
} else if (case=="IV") {
 Z1<- rbind(Z1,1:T)
 rownames(Z1)[m+1]<- "t"
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
if (0) {V <- orthonormalization(V,basis=F,norm=T)}

## compute test statistics for all ranks 0:n ########################
## likelihood ratio test trace statistic: H(r) in H(n)
# -2*log Q(H(r)|H(n))
LR.trace<- vector("list",n); names(LR.trace)<- paste("rank",0:(n-1),"vs.",n)
for (r in 0:(n-1) ) {
 LR.trace[[paste("rank",r,"vs.",n)]]<- -T*sum(log(1-lambda[(r+1):n]))
}
## likelihood ratio test max eigenvalue statistic: H(r) in H(r+1)
# -2*log Q(H(r)|H(r+1))
LR.maxeigen<- vector("list",n); names(LR.maxeigen)<- paste("rank",0:(n-1),"vs.",1:n)
for (r in 0:(n-1) ) {
 LR.maxeigen[[paste("rank",r,"vs.",r+1)]]<- -T*log(1-lambda[r+1])
}
## give critical values of trace statistic for all ranks 0:n
# critical values are taken from [1]
CV.trace<- vector("list",n); names(CV.trace)<- paste("rank",0:(n-1),"vs.",n)
for (r in 0:(n-1) ) {
  cv5<-  CV.trace.table[[paste("case",case,"5%")]][13-(n-r),k+1]
  #cv10<- CV.trace.table[[paste("case",case,"10%")]][13-(n-r),k+1]
  #cv<- matrix(c(cv10,cv5),ncol=2,nrow=1)
  #colnames(cv)<-c("10%","5%")
  cv<- matrix(cv5,ncol=1,nrow=1)
  colnames(cv)<-c("5%")
  CV.trace[[paste("rank",r,"vs.",n)]]<-cv
 }
CV.maxeigen<- vector("list",n); names(CV.maxeigen)<- paste("rank",0:(n-1),"vs.",1:n)
for (r in 0:(n-1) ) {
  cv5<-  CV.maxeigen.table[[paste("case",case,"5%")]][13-(n-r),k+1]
  #cv10<- CV.maxeigen.table[[paste("case",case,"10%")]][13-(n-r),k+1]
  #cv<- matrix(c(cv10,cv5),ncol=2,nrow=1)
  #colnames(cv)<-c("10%","5%")
  cv<- matrix(cv5,ncol=1,nrow=1)
  colnames(cv)<-c("5%")
  
  CV.maxeigen[[paste("rank",r,"vs.",r+1)]]<-cv
 }

## output:
mdls<-list(type="weakly exogenous VECM",freq=freq,m=m,n=n,p=p,q=q,T=T,case=case,lambda=lambda,LR.trace=LR.trace,CV.trace=CV.trace,LR.maxeigen=LR.maxeigen,CV.maxeigen=CV.maxeigen)
return(mdls)
}

