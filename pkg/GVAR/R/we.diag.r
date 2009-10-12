we.diag <- function(Z,n,p,q=p,case)
   #     Z ... matrix of endogenous (Y) and weakly exogenous (X) time series
   #     n ... number of Y variables
   #     p ... lag order of endogenous variables
   #     q ... lag order of weakly exogenous variables
   #  case ... case I-V with respect to PSS 2000
{

## --- Diagnostic tests for the weakly exogenous varibales X ---
#
#  The exogenous variables X must be (a) integrated of order 1 (X~I(1));
#  it follows (b) that X is not cointegrated on its own, and (c) that the
#  differenced process does not depend on the lagged Z.

  source("cv.tables.pesaran.r")
  
  m <- dim(Z)[2]                        # number of total variables
  k <- m-n                              # number of X variables
  T.data <- dim(Z)[1]
  T <- T.data-max(p,q)
  
  X <- Z[,(n+1):m]
  Y <- Z[,1:n]
  DZ <- diff(Z)
  DY <- DZ[max(p,q):(T.data-1),1:n]
  DX <- DZ[max(p,q):(T.data-1),(n+1):(n+k)]

  if (max(p,q)>1)
  {
    DZstar <- NULL
    if (min(p,q)>1)
    {
      for (i in 1:(min(p,q)-1))
      {
        DZstar <- cbind(DZstar,DZ[(max(p,q)-i):(T.data-(i+1)),])
      }
    }
    if (p>q)
    {
      for (i in q:(p-1))
      {
        DZstar <- cbind(DZstar,DZ[(max(p,q)-i):(T.data-(i+1)),1:n])
      }
    } else if (q>p) {
      for (i in p:(q-1))
      {
        DZstar <- cbind(DZstar,DZ[(max(p,q)-i):(T.data-(i+1)),(n+1):(n+k)])
      }
    }
  } else {
    DZstar <- rep(0,T)
  }

  Ones <- rep(1,T)
  Trend <- 1:T

# --- Generate Z* and Design matrix for the different cases

  if (case=="I")
  {
    Zstar <- Z[max(p,q):(T.data-1),]
    Design_woX <- DZstar
    lNewn <- n
  }
  if (case=="II")
  {
    Zstar <- cbind(Z[max(p,q):(T.data-1),],Ones)
    Design_woX <- DZstar
    lNewn <- n+1
  }
  if (case=="III")
  {
    Zstar <- Z[max(p,q):(T.data-1),]
    Design_woX <- cbind(Ones,DZstar)
    lNewn <- n
  }
  if (case=="IV")
  {
    Zstar <- cbind(Z[max(p,q):(T.data-1),],Trend)
    Design_woX <- cbind(Ones,DZstar)
    lNewn <- n+1
  }
  if (case=="V")
  {
    Zstar <- Z[max(p,q):(T.data-1),]
    Design_woX <- cbind(Ones,Trend,DZstar)
    lNewn <- n
  }


  Design <- cbind(DX,Design_woX)

  if (max(p,q)==1){Design <- Design[,1:(dim(Design)[2]-1)]}

# --- Compute projector matrix with singular value decomposition (SVD) to compute residuals

  SVD <- svd(Design)
  V <- SVD$v
  U <- SVD$u
  S <- diag(SVD$d)
  lNbCols <- dim(V)[2]
  Uproject <- U[,1:lNbCols]
  Projector <- tcrossprod(Uproject)

  Zhat <- Zstar - Projector%*%Zstar
  DYhat <- DY - Projector%*%DY

  Syy <- (t(DYhat)%*%DYhat)/T
  Syz <- (t(DYhat)%*%Zhat)/T
  Szz <- (t(Zhat)%*%Zhat)/T
  Sall <- cbind(rbind(Syy,t(Syz)),rbind(Syz,Szz))

  m_ <- lNewn+k
  L <- t(chol(Sall))
  Lz <- t(chol(Sall[(n+1):(m_+n),(n+1):(m_+n)]))
  L21 <- L[(n+1):(m_+n),1:n]
  LzInv <- solve(Lz)

# --- Compute Eigenvalues/-vectors

  eig <- eigen(LzInv%*%L21%*%t(L21)%*%t(LzInv))
  EigvecsFinalTranspo <- LeftEigvecs <- t(LzInv)%*%eig$vectors
  EigvalsFinal <- eig$values

## --- Exogeneity diagnostics ---

# --- Test if DX depends on Y(t-1)

  CointResid <- Zstar%*%t(EigvecsFinalTranspo)
  if (max(p,q)>1 || (case!="I" && case!= "II"))
  {
    if (max(p,q)==1){Xdesign <- Design_woX[,1:(dim(Design_woX)[2]-1)]}else{Xdesign <- Design_woX}
    SVD2 <- svd(Xdesign)
    V <- SVD2$v
    U <- SVD2$u
    lNbCols <- dim(V)[2]
    Uproject <- U[,1:lNbCols]
    Projector <- Uproject%*%t(Uproject)
  }else
  {
    Projector <- matrix(0,nrow=T,ncol=T)
  }

  ExRestr <- DX-Projector%*%DX
  CointResidhat <- CointResid-Projector%*%CointResid

  alpha_xy <- list()
  Test_alpha_xy <- vector()
  critval.alpha.xy <- vector()
  prob.alpha.xy <- vector()

  for (r_ in 1:n)
  {
    SVDc <- svd(CointResidhat[,1:r_])
    Vc <- SVDc$v
    Uc <- SVDc$u
    Sc <- SVDc$d
    if (r_ == 1)
    {
      alpha_xy[[r_]] <- Vc%*%(1/Sc)%*%t(Uc[,1:r_])%*%ExRestr
    }
    if (r_ > 1)
    {
      alpha_xy[[r_]] <- Vc%*%diag(1/Sc)%*%t(Uc[,1:r_])%*%ExRestr
    }
    ExFree <- ExRestr-CointResidhat[,1:r_]%*%alpha_xy[[r_]]
    Test_alpha_xy[r_] <- T*(log(det(t(ExRestr)%*%ExRestr))-log(det(t(ExFree)%*%ExFree)))
    critval.alpha.xy[r_] <- qchisq(.95, df=k*r_)
    prob.alpha.xy[r_] <- pchisq(Test_alpha_xy[r_],df=k*r_,lower.tail = F)       # Null: X independent of Y levels
  }

# --- Test if X is cointegrated

  Xstar <- Zstar[,(n+1):m_]
  DXhat <- ExRestr
  if (q>1){XhatStar <- Xstar-Projector%*%Xstar}else{XhatStar <- Xstar}

  S00 <- t(DXhat)%*%DXhat/T
  S01 <- t(DXhat)%*%XhatStar/T
  S11 <- t(XhatStar)%*%XhatStar/T
  SXcoint <- cbind(rbind(S00,t(S01)),rbind(S01,S11))

  m_ <- max(dim(SXcoint))
  n_ <- k

  L.x <- t(chol(SXcoint))
  Lz.x <- t(chol(SXcoint[(n_+1):m_,(n_+1):m_]))
  L21.x <- L.x[(n_+1):m_,1:n_]
  LzInv.x <- solve(Lz.x)

  eig.x <- eigen(LzInv.x%*%L21.x%*%t(L21.x)%*%t(LzInv.x))
  EigvecsXcointTranspo <- LeftEigvecs.x <- t(LzInv.x)%*%eig.x$vectors
  EigvalsXcoint <- eig.x$values

  LogEigsX <- -T*log(rep(1,length(EigvalsXcoint))-EigvalsXcoint)

# --- Print output

  TestStata <- cbind(round(sum(LogEigsX),4),CV.trace.table[[paste("case",case,"5%")]][13-k,1])
  rownames(TestStata) <- "X ~ I(1)           "
  colnames(TestStata) <- c("    teststat","    crit.value")
  TestStatb <- cbind(round(LogEigsX[1],4),CV.maxeigen.table[[paste("case",case,"5%")]][13-k,1])
  rownames(TestStatb) <- "No Cointegration   "
  colnames(TestStatb) <- c("    teststat","    crit.value")
  TestStatc <- cbind(round(Test_alpha_xy,4),round(prob.alpha.xy,4))
  rownames(TestStatc) <- paste("r = ",1:n,":             ",sep="")
  colnames (TestStatc) <- c("    teststat","       p-value")

  cat("\n------------ Exogeneity Diagnostics ------------\n\n")
  cat("a. Test if X is stationary:\n")
  print(TestStata)
  cat("\n\n")
  cat("b. Test if X is cointegrated:\n")
  print(TestStatb)
  cat("\n\n")
  cat("c. Test if lagged cointegration relationship affects DX:\n")
  print(TestStatc)
  cat("-------------------------------------------------\n")
}


