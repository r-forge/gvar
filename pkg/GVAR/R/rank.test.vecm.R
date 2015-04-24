rank.test.vecm <- function (Y.ts,etw,p,case,ex=0,lex=NULL,season=NULL,season.start.time=NULL)
## notation: almost as in Johansen's book "Likelihood-based Inference in Cointegrated VAR Models"
## model:
# \Delta Y_t = \Pi Y_{t-1}+\sum_{i=1}^{k-1} \Gamma_i\Delta Y_{t-i}+\mu_0+\mu_1 t+\Phi D_t+\epsilon_t, where \epsilon_t is N(0,\Omega)

{
freq<- etw[["freq"]] # time sampling frequency
dt<- 1/freq # time sampling interval
T<- (etw[["end"]]-etw[["start"]])*freq+1  # number of time samples for estimation
n<- dim(Y.ts)[2]-ex # number of variables in Y_t

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
 rownames(Dt)<- c("constant","t")
} else if (case=="H_1(r)") {
 Dt<- rbind(rep(1,T)); rownames(Dt)<- c("constant")
} else if (case=="H_2(r)") {
 Dt<- rbind();
} else if (case=="H^*(r)") {
 Dt<- rbind(rep(1,T),seq(etw[["start"]],etw[["end"]],by=dt))
 rownames(Dt)<- c("constant","t")
} else if (case=="H_1^*(r)") {
 Dt<- rbind(rep(1,T)); rownames(Dt)<- c("constant")
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
 rownames(Z2)[(p-1)*n+1]<- "constant"
} else if (case=="H_1^*(r)") {
 Z1<- rbind(Z1, 1)
 rownames(Z1)[n+1]<- "constant"
} else {stop("\nUnkown case.\n")}
if (!(is.null(season))) { # seasonal dummies
 l<- (season.start.time*freq-etw[["start"]]*freq)%%season
 dum<- (diag(season) - 1/season)[-season,]
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

## compute test statistics for all ranks #############################
## likelihood ratio test trace statistic: H(r) in H(n) for all ranks
# -2*log Q(H(r)|H(n))
LR.trace<- vector("list", n); names(LR.trace)<- paste("rank",0:(n-1),"vs.",n)
for (r in 0:(n-1) ) {
 LR.trace[[paste("rank",r,"vs.",n)]]<- -T*sum(log(1-lambda[(r+1):n]))
}

## likelihood ratio test max eigenvalue statistic: H(r) in H(r+1) for all ranks
# -2*log Q(H(r)|H(r+1))
LR.maxeigen<- vector("list", n); names(LR.maxeigen)<- paste("rank",0:(n-1),"vs.",1:n)
for (r in 0:(n-1) ) {
 LR.maxeigen[[paste("rank",r,"vs.",r+1)]]<- -T*log(1-lambda[r+1])
}

## give critical values of trace statistic for all ranks #############
# critical values are taken from Johansen's book "Likelihood-based Inference in Cointegrated VAR Models"
if (1) {
 data(cv.tables.johansen)
 if (case=="H(r)") {          # Table 15.5
  CV.trace.table <- cv.tables.johansen$table15.5
 } else if (case=="H_1(r)") {   # Table 15.3
  CV.trace.table <- cv.tables.johansen$table15.3
 } else if (case=="H_2(r)") {   # Table 15.1
  CV.trace.table <- cv.tables.johansen$table15.1
 } else if (case=="H^*(r)") {   # Table 15.4
  CV.trace.table <- cv.tables.johansen$table15.4
 } else if (case=="H_1^*(r)") { # Table 15.2
  CV.trace.table <- cv.tables.johansen$table15.2
 }
 CV.trace<- vector("list", n); names(CV.trace)<- paste("rank",0:(n-1),"vs.",n)
 for (r in 0:(n-1) ) {
  p_r<- n-r
  CV.trace[[paste("rank",r,"vs.",n)]]<- CV.trace.table[p_r,]
 }
}

## output ##############################################################
# estimates of: alpha, beta, Pi, Gamma, mu_t, Omega
# t-values of beta, cf. p.100 in Luetkepohl
# various test statistics and critical values

# browser()
# collect all estimated models for r=0:n in a list:
mdls <- list(type="pure VECM",freq=freq,n=n,p=p,T=T,season=season,season.start.time=season.start.time,case=case,lambda=lambda,LR.trace=LR.trace,CV.trace=CV.trace,LR.maxeigen=LR.maxeigen)
return(mdls)
}
