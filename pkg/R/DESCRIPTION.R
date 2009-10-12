summary.vecm <- function(object, ...)
{
  if (object[["type"]]=="pure VECM") 
  {
    roundto <- 5
    alpha <- round(object$alpha,roundto)
    beta <- round(object$beta,roundto)
    rownames(beta) <- colnames(object$dat)
    colnames(beta) <- paste(1:object$r,".:   ",sep="")
    Beta <- NULL
    Gamma <- object$Gamma
    const <- cbind(object$mu0,object$mu1)
    nam <- vector()
    if (!is.null(object$mu0)) {nam <- c(nam,"         Const")}
    if (!is.null(object$mu1)) {nam <- c(nam,"Trend")}
    if (!is.null(const))
    {
      const <- round(const,roundto)
      colnames(const) <- nam
    }
    if (!is.null(object$season)) 
    {
      Phi <- round(object$Phi,roundto)
      nam <- "            s1"
      if (object$season>1) {nam <- c(nam,paste("s",2:(object$season-1),sep=""))}
      colnames(Phi) <- nam 
    }
  
    cat("\nTime series information:\n")
    cat("Cointegration Rank: ",object$r,"\t\tLag Order: ",object$p,"\n")
    cat("\n")
    cat("Coefficients:\n")

    cat("Beta':\n")
    tvalsB <- round(object$tvals$beta,rounto)
    l <- 0
    for (i in 1:object$r)
    {
      temptB <- rep("-",object$r)
      seTemp <- rep("-",object$r)
      for (j in l+(1:(object$n-object$r)))
      {
        temptB <- c(temptB,paste("[",round(tvalsB[j],2),"]",sep=""))        
      }
      seTemp <- c(seTemp,rep("NA",(object$n-object$r)))

      BETA <- rbind(beta[,i],seTemp,temptB)
      rownames(BETA) <- c(colnames(beta)[i],"(Std.Err.)","[t-Value]")
      print(as.data.frame(BETA))
      l <- i*(object$n-object$r)
    }
    cat("\n")

    cat("Alpha:\n")
    print(alpha)
    cat("\n")
    if (object$p>1)
    {
      Gammer <- NULL
      seG <- NULL
      tvalsG <- NULL
      pvalsG <- NULL
      cat("Gamma:\n")
      for (i in 1:length(object$Gamma))
      {        
        colnames(Gamma[[i]])[1] <- paste("        ",colnames(Gamma[[i]])[1])
        Gammer <- cbind(Gammer,round(Gamma[[i]],roundto))
        seG <- cbind(seG,round(object$se$Gamma[[i]],roundto))
        tvalsG <- cbind(tvalsG,round(object$tvals$Gamma[[i]],roundto))
        pvlasG <- cbind(pvalsG,round(object$pvals$Gamma[[i]],roundto))
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:length(Gammer[i,]))
        {
          seTemp <- c(seTemp,paste("(",round(seG[i,j],3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",round(tvalsG[i,j],2),"]",sep=""))
        }
        GAMMA <- rbind(Gammer[i,],seTemp,tvalTemp)
        rownames(GAMMA) <- c(rownames(Gammer)[i],"(Std.Err.)","[t-Value]")
        print(as.data.frame(GAMMA))
      }
      cat("\n")
    }
    if (!is.null(object$mu0) || !is.null(object$mu1))
    { 
      cat("Intercept (and Trend) in VAR:\n")
      if (object$case=="V")
      {
        seC <- cbind(object$se$mu0,object$se$mu1)
        tvalC <- cbind(object$se$mu0,object$se$mu1)
        pvalC <- cbind(object$se$mu0,object$se$mu1)
      } else if (object$case=="IV" || object$case=="H_1(r)"){
        seC <- cbind(object$se$mu0)
        tvalC <- cbind(object$se$mu0)
        pvalC <- cbind(object$se$mu0)
      } else {
        seC <- NULL
        tvalC <- NULL
        pvalC <- NULL
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:dim(const)[2])
        {
          seTemp <- c(seTemp,paste("(",round(seC[i,j],3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",round(tvalC[i,j],2),"]",sep=""))
        }
        CONST <- rbind(const[i,],seTemp,tvalTemp)
        rownames(CONST) <- c(rownames(const)[i],"(Std.Err.)","[t-Value]")
        colnames(CONST) <- colnames(const)
        print(as.data.frame(CONST))
      }
      cat("\n")                            
    }

    if (!is.null(object$season))
    {
      cat("Seasonality:\n")
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:(object$season-1))
        {
          seTemp <- c(seTemp,paste("(",round(object$se$season[i,j],3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",round(object$tvals$season[i,j],2),"]",sep=""))
        }
        SEASON <- rbind(Phi[i,],seTemp,tvalTemp)
        rownames(SEASON) <- c(rownames(object$se$season)[i],"(Std.Err.)","[t-Value]")
        colnames(SEASON) <- colnames(Phi)
        print(as.data.frame(SEASON))
      }
      cat("\n")
    }
    cat("\n")
  } else if (object[["type"]]=="weakly exogenous VECM") {

  }
}