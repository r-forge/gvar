summary.vecm <- function(object, ...)
{
  if (object[["type"]]=="pure VECM") 
  {
    roundto <- 5
    if (!(object$r==object$n))
    {
      alpha <- round(object$alpha,roundto)
      beta <- round(object$beta,roundto)
      colnames(beta) <- paste(1:object$r,".:   ",sep="")
      Beta <- NULL
    }
    Gamma <- object$Gamma
    Psi <- object$Psi
    const <- cbind(object$mu0,object$mu1)
    nam <- vector()
    if (!is.null(object$mu0)) {nam <- c(nam,"       Const")}
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
	if (!(object$r==object$n))
	{ 
      cat("Beta':\n")
      if (object$r==1)
      {
        tvalsB <- matrix(round(object$tvals$beta,roundto),ncol=1)
      } else {
        tvalsB <- round(object$tvals$beta,roundto)
      }
      seB <- round(object$se$beta,roundto)
      for (i in 1:object$r)
      {
        temptB <- rep("-",object$r)
        seTemp <- rep("-",object$r)        
        if (object$case=="H_1^*(r)" || object$case=="H^*(r)")
        {
      	  i.b <- object$n-object$r+1
        } else {
      	  i.b <- object$n-object$r
        }
        for (j in 1:i.b)
        {
          seTemp <- c(seTemp,paste("(",format(round(seB[j,i],3),nsmall=3),")",sep=""))
          temptB <- c(temptB,paste("[",format(round(tvalsB[j,i],2),nsmall=2),"]",sep=""))        
        }

        BETA <- rbind(format(beta[,i],nsmall=roundto),seTemp,temptB)
        BETA <- rbind(BETA,rep("",dim(BETA)[2]))
        rownames(BETA) <- c(colnames(beta)[i],"(Std.Err.)","[t-Value]","")
        print(as.data.frame(BETA))
      }
      cat("\n")        

      cat("Alpha:\n")
      print(alpha)
      cat("\n")
    }
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
          seTemp <- c(seTemp,paste("(",format(round(seG[i,j],3),nsmall=3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",format(round(tvalsG[i,j],2),nsmall=2),"]",sep=""))
        }
        GAMMA <- rbind(format(Gammer[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(GAMMA) <- c(rownames(Gammer)[i],"(Std.Err.)","[t-Value]")
        print(as.data.frame(GAMMA))
      }
      cat("\n")
    }
    if (object$ex>0)
    {
      Psis <- NULL
      sePs <- NULL
      tvalsPs <- NULL
      pvalsPs <- NULL
      cat("Psi:\n")
      for (i in 1:length(object$Psi))
      {        
        colnames(Psi[[i]])[1] <- paste("        ",colnames(Psi[[i]])[1])
        Psis <- cbind(Psis,round(Psi[[i]],roundto))
        sePs <- cbind(sePs,round(object$se$Psi[[i]],roundto))
        tvalsPs <- cbind(tvalsPs,round(object$tvals$Psi[[i]],roundto))
        pvlasPs <- cbind(pvalsPs,round(object$pvals$Psi[[i]],roundto))
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:length(Psis[i,]))
        {
          seTemp <- c(seTemp,paste("(",format(round(sePs[i,j],3),nsmall=3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",format(round(tvalsPs[i,j],2),nsmall=2),"]",sep=""))
        }
        PSI <- rbind(format(Psis[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(PSI) <- c(rownames(Psis)[i],"(Std.Err.)","[t-Value]")
        print(as.data.frame(PSI))
      }
      cat("\n")
    }
    if (!is.null(object$mu0) || !is.null(object$mu1))
    { 
      cat("Intercept (and Trend) in VECM:\n")
      if (object$case=="H(r)")
      {
        seC <- cbind(object$se$mu0,object$se$mu1)
        tvalC <- cbind(object$tvals$mu0,object$tvals$mu1)
        pvalC <- cbind(object$pvals$mu0,object$pvals$mu1)
        n.case <- 2
      } else if (object$case=="H^*(r)" || object$case=="H_1(r)"){
        seC <- cbind(object$se$mu0)
        tvalC <- cbind(object$tvals$mu0)
        pvalC <- cbind(object$pvals$mu0)
        n.case <- 1
      } else {
        seC <- NULL
        tvalC <- NULL
        pvalC <- NULL
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        
        if (object$case!="H_1^*(r)")
        {
          for (j in 1:n.case)
          {
            seTemp <- c(seTemp,paste("(",format(round(seC[i,j],3),nsmall=3),")",sep=""))
            tvalTemp <- c(tvalTemp,paste("[",format(round(tvalC[i,j],2),nsmall=2),"]",sep=""))
          }
        }
        if (object$case=="H_1^*(r)" || object$case=="H^*(r)")
        {
          seTemp <- c(seTemp,"-")
          tvalTemp <- c(tvalTemp,"-")
        } 
        
        CONST <- rbind(format(const[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(CONST) <- c(rownames(const)[i],"(Std.Err.)","[t-Value]")
        colnames(CONST) <- paste("   ",colnames(const),sep="")
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
          seTemp <- c(seTemp,paste("(",format(round(object$se$season[i,j],3),nsmall=3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",format(round(object$tvals$season[i,j],2),nsmall=2),"]",sep=""))
        }
        SEASON <- rbind(format(Phi[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(SEASON) <- c(rownames(object$se$season)[i],"(Std.Err.)","[t-Value]")
        colnames(SEASON) <- paste("   ",colnames(Phi),sep="")
        print(as.data.frame(SEASON))
      }
      cat("\n")
    }
    cat("\n")
  } else if (object[["type"]]=="weakly exogenous VECM") {
  	
  	
  	#--------------- WEAKLY EXOGENOUS VECM --------------
  	
    roundto <- 5
    if (!(object$r==object$n))
    {
      alpha <- round(object$alpha,roundto)
      beta <- round(object$beta,roundto)
      if (object$case=="II") 
      {
         rownames(beta) <- c(colnames(object$dat[,1:object$m]),"Const")
      } else if (object$case=="IV") {
         rownames(beta) <- c(colnames(object$dat[,1:object$m]),"Trend")
      } else {     
    	    rownames(beta) <- colnames(object$dat[,1:object$m])
      }
      colnames(beta) <- paste(1:object$r,".:   ",sep="")
      Beta <- NULL
    }
    Lambda <- object$Lambda
    Phi <- object$Phi
    Psi <- object$Psi
    const <- cbind(object$c.0,object$c.1)
    nam <- vector()
    if (!is.null(object$c.0)) {nam <- c(nam,"       Const")}
    if (!is.null(object$c.1)) {nam <- c(nam,"Trend")}
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
    cat("Cointegration Rank: ",object$r,
        "\nLag Order of endogenous Variables: ",object$p,
        "\nLag Order of weakly exogenous Variables: ",object$q,"\n")
    cat("\n")
    cat("Coefficients:\n")

	if (!(object$r==object$n))
	{
      cat("Beta':\n")
      if (object$r==1)
      {
        tvalsB <- matrix(round(object$tvals$beta,roundto),ncol=1)
      } else {
        tvalsB <- round(object$tvals$beta,roundto)
      }
      seB <- round(object$se$beta,roundto)
      for (i in 1:object$r)
      {
        temptB <- rep("-",object$r)
        seTemp <- rep("-",object$r)        
        if (object$case=="II" || object$case=="IV")
        {
      	  i.b <- object$m-object$ex-object$r+1
        } else {
      	  i.b <- object$m-object$ex-object$r
        }
        for (j in 1:i.b)
        {
          seTemp <- c(seTemp,paste("(",format(round(seB[j,i],3),nsmall=3),")",sep=""))
          temptB <- c(temptB,paste("[",format(round(tvalsB[j,i],2),nsmall=2),"]",sep=""))        
        }

        BETA <- rbind(format(beta[,i],nsmall=roundto),seTemp,temptB)
        rownames(BETA) <- c(colnames(beta)[i],"(Std.Err.)","[t-Value]")
        print(as.data.frame(BETA))
        cat("\n")
      }
#      cat("\n")

      cat("Alpha:\n") 
      colnames(alpha)[1] <- paste("        ",colnames(alpha)[1])
      print(alpha)
      cat("\n")
    }
    cat("Lambda:\n")
    colnames(Lambda)[1] <- paste("        ",colnames(Lambda)[1])
    seL <- round(object$se$Lambda,roundto)
    tvalsL <- round(object$tvals$Lambda,roundto)
    pvalsL <- round(object$pvals$Lambda,roundto)
    for (i in 1:object$n)
    {
      seTemp <- NULL
      tvalTemp <- NULL
      for (j in 1:(dim(object$dat)[2]-object$n))
      {
        seTemp <- c(seTemp,paste("(",format(round(seL[i,j],3),nsmall=3),")",sep=""))
        tvalTemp <- c(tvalTemp,paste("[",format(round(tvalsL[i,j],2),nsmall=2),"]",sep=""))
      }
      LAMBDA <- rbind(format(round(Lambda[i,],roundto),nsmall=roundto),seTemp,tvalTemp)
      rownames(LAMBDA) <- c(rownames(Lambda)[i],"(Std.Err.)","[t-Value]")
      if (dim(Lambda)[2]<2){colnames(LAMBDA) <- colnames(Lambda)}
      print(as.data.frame(LAMBDA))
    }
    cat("\n")
    if (object$p>1)
    {
      tempPhi <- NULL
      sePh <- NULL
      tvalsPh <- NULL
      pvalsPh <- NULL
      cat("Phi:\n")
      for (i in 1:length(object$Phi))
      {        
        colnames(Phi[[i]])[1] <- paste("        ",colnames(Phi[[i]])[1])
        tempPhi <- cbind(tempPhi,round(Phi[[i]],roundto))
        sePh <- cbind(sePh,round(object$se$Phi[[i]],roundto))
        tvalsPh <- cbind(tvalsPh,round(object$tvals$Phi[[i]],roundto))
        pvlasPh <- cbind(pvalsPh,round(object$pvals$Phi[[i]],roundto))
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:length(tempPhi[i,]))
        {
          seTemp <- c(seTemp,paste("(",format(round(sePh[i,j],3),nsmall=3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",format(round(tvalsPh[i,j],2),nsmall=2),"]",sep=""))
        }
        PHI <- rbind(format(tempPhi[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(PHI) <- c(rownames(tempPhi)[i],"(Std.Err.)","[t-Value]")
        print(as.data.frame(PHI))
      }
      cat("\n")
    }
    if ((object$q>1) || (object$lex>1))
    {
      tempPsi <- NULL
      sePs <- NULL
      tvalsPs <- NULL
      pvalsPs <- NULL
      cat("Psi:\n")
      for (i in 1:length(object$Psi))
      {        
        colnames(Psi[[i]])[1] <- paste("        ",colnames(Psi[[i]])[1])
        tempPsi <- cbind(tempPsi,round(Psi[[i]],roundto))
        sePs <- cbind(sePs,round(object$se$Psi[[i]],roundto))
        tvalsPs <- cbind(tvalsPs,round(object$tvals$Psi[[i]],roundto))
        pvlasPs <- cbind(pvalsPs,round(object$pvals$Psi[[i]],roundto))
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        for (j in 1:length(tempPsi[i,]))
        {
          seTemp <- c(seTemp,paste("(",format(round(sePs[i,j],3),nsmall=3),")",sep=""))
          tvalTemp <- c(tvalTemp,paste("[",format(round(tvalsPs[i,j],2),nsmall=2),"]",sep=""))
        }
        PSI <- rbind(format(tempPsi[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(PSI) <- c(rownames(tempPsi)[i],"(Std.Err.)","[t-Value]")
        if (dim(Psi[[1]])[2]<2){colnames(PSI) <- colnames(Psi[[1]])}
        print(as.data.frame(PSI))
      }
      cat("\n")
    }

    if (!is.null(object$c.0) || !is.null(object$c.1))
    { 
      cat("Intercept (and Trend) in VECM:\n")
      if (object$case=="V")
      {
        seC <- cbind(object$se$c.0,object$se$c.1)
        tvalC <- cbind(object$tvals$c.0,object$tvals$c.1)
        pvalC <- cbind(object$pvals$c.0,object$pvals$c.1)
        n.case <- 2
      } else if (object$case=="IV" || object$case=="III"){
        seC <- cbind(object$se$c.0)
        tvalC <- cbind(object$tvals$c.0)
        pvalC <- cbind(object$pvals$c.0)
        n.case <- 1
      } else {
        seC <- NULL
        tvalC <- NULL
        pvalC <- NULL
      }
      for (i in 1:object$n)
      {
        seTemp <- NULL
        tvalTemp <- NULL
        if (object$case!="II")
        {
          for (j in 1:n.case)
          {
            seTemp <- c(seTemp,paste("(",format(round(seC[i,j],3),nsmall=3),")",sep=""))
            tvalTemp <- c(tvalTemp,paste("[",format(round(tvalC[i,j],2),nsmall=2),"]",sep=""))
          }
        }
        if (object$case=="II" || object$case=="IV")
        {
          seTemp <- c(seTemp,"-")
          tvalTemp <- c(tvalTemp,"-")
        } 
        CONST <- rbind(format(const[i,],nsmall=roundto),seTemp,tvalTemp)
        rownames(CONST) <- c(rownames(const)[i],"(Std.Err.)","[t-Value]")
        colnames(CONST) <- paste("    ",colnames(const),sep="")
        print(as.data.frame(CONST))
      }
      cat("\n")                            
    }
  }
}