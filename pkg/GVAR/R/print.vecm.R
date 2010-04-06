print.vecm <- function(x, ...)
{
if (x[["type"]]=="pure VECM") 
  {
    roundto <- 5
    if (!(x$r==x$n))
    {
      alpha <- round(x$alpha,roundto)
      beta <- round(x$beta,roundto)
      colnames(beta) <- paste(1:x$r,".:   ",sep="")
    }
    Gamma <- x$Gamma
    const <- cbind(x$mu0,x$mu1)
    nam <- vector()
    if (!is.null(x$mu0)) {nam <- c(nam,"         Const")}
    if (!is.null(x$mu1)) {nam <- c(nam,"Trend")}
    if (!is.null(const))
    {
      const <- round(const,roundto)
      colnames(const) <- nam
    }
    if (!is.null(x$season)) 
    {
      Phi <- round(x$Phi,roundto)
      nam <- "            s1"
      if (x$season>1) {nam <- c(nam,paste("s",2:(x$season-1),sep=""))}
      colnames(Phi) <- nam 
    }  
  
    cat("\n")
    cat("Coefficients:\n")
    if (!(x$r==x$n))
    {
      cat("Beta:\n")
      print(t(beta))
      cat("\n")
      cat("Alpha:\n")
      print(alpha)
      cat("\n")
    }
    if (x$p>1)
    {
      Gammer <- NULL
      cat("Gamma:\n")
      for (i in 1:length(x$Gamma))
      {        
        colnames(Gamma[[i]])[1] <- paste("        ",colnames(Gamma[[i]])[1])
        Gammer <- cbind(Gammer,round(Gamma[[i]],roundto))
      }
      print(Gammer)
      cat("\n")
    }
    if (!is.null(x$mu0) || !is.null(x$mu1))
    { 
      cat("Intercept (and Trend) in VAR:\n")
      print(const)
      cat("\n")
    }

    if (!is.null(x$season))
    {
      cat("Seasonality:\n")
      print(Phi)
      cat("\n")
    }
    cat("\n")
    
    
    
  } else if (x[["type"]]=="weakly exogenous VECM") {
    roundto <- 5
    if (!(x$r==x$n))
    {
      alpha <- round(x$alpha,roundto)
      beta <- round(x$beta,roundto)
      if (x$case=="II") 
      {
         rownames(beta) <- c(colnames(x$dat[,1:x$m]),"Const")
      } else if (x$case=="IV") {
         rownames(beta) <- c(colnames(x$dat[,1:x$m]),"Trend")
      } else {     
    	    rownames(beta) <- colnames(x$dat[,1:x$m])
      }
      colnames(beta) <- paste(1:x$r,".:   ",sep="")
    }
    Lambda <- round(x$Lambda,roundto)
    const <- cbind(x$c.0,x$c.1)
    nam <- vector()
    if (!is.null(x$c.0)) {nam <- c(nam,"         Const")}
    if (!is.null(x$c.1)) {nam <- c(nam,"Trend")}
    if (!is.null(const))
    {
      const <- round(const,roundto)
      colnames(const) <- nam
    }
    if (!is.null(x$Phi)){Phi <- x$Phi}
    if (!is.null(x$Psi)){Psi <- x$Psi}  
  
    cat("\n")
    cat("Coefficients:\n")
    if (!(x$r==x$n))
    {
      cat("Beta:\n")
      print(t(beta))
      cat("\n")
      cat("Alpha:\n")
      colnames(alpha)[1] <- paste("        ",colnames(alpha)[1])
      print(alpha)
      cat("\n")    
    }
    cat("Lambda:\n")
    colnames(Lambda)[1] <- paste("        ",colnames(Lambda)[1])
    print(Lambda)
    cat("\n")
    if (x$p>1)
    {
      PHI <- NULL
      cat("Phi:\n")
      for (i in 1:length(x$Phi))
      {        
        colnames(Phi[[i]])[1] <- paste("        ",colnames(Phi[[i]])[1])
        PHI <- cbind(PHI,round(Phi[[i]],roundto))
      }
      print(PHI)
      cat("\n")
    }
    if ((x$q>1) || (x$lex>1))
    {
      PSI <- NULL
      cat("Psi:\n")
      for (i in 1:length(x$Psi))
      {        
        colnames(Psi[[i]])[1] <- paste("        ",colnames(Psi[[i]])[1])
        PSI <- cbind(PSI,round(Psi[[i]],roundto))
      }
      print(PSI)
      cat("\n")
    }
    if (!is.null(x$c.0) || !is.null(x$c.1))
    { 
      cat("Intercept (and Trend) in VAR:\n")
      print(const)
      cat("\n")
    }
#
#    if (!is.null(x$season))
#    {
#      cat("Seasonality:\n")
#      print(Phi)
#      cat("\n")
#    }
    cat("\n")
  }
}