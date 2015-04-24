# Forecast Error Variance Decomposition as in Peseran and Shin 1996
#
# input
# x         ...  a gvar object
# n         ...  number of periods, in the GVAR literature often set to 40
# shock.var ...  vector of length 2, first element specifies the country to be shocked, second element the variable in the country (e.g. shock.var=c(13,3) third variable in US is going to be shocked)
# shock.dir ...  a character, if set to "-1", the shock is going to be a negative 1sd shock, where sd refers to the country specific residuals. Alternatively, you can submit any number indicating
#                that the shock  has to be standardized, for example a variable that is measured in % might get shocked by "-0.03/sd"
#
# output
# psi       ...  matrix of impulse responses
# Fmat      ...  coefficient matrix of ar1 representation of the system, important to check stability of system via its eigenvalues
# G         ...  G
# H         ...  H
# U         ...  U (residual matrix of the gvar system)
# shock.dir ...  contains the submitted shock (either in standard deviations or in other magnitudes that get standardized within GIR)
# sigma.il  ...  residual standard deviation (i.e. magnitude of the shock if not elsewise specified)

FEVD=function (x, n=40, shock.var, shock.dir=-1,scal=FALSE){

    endoN=sapply(x$we.vecms,function(x) x$n)
    G <- x$G
    U <- x$U
    H.list <- x$H
    H <- H.list[[1]]
    p = x$arguments$p
    q = x$arguments$q
    l <- length(x$subsys)
    if (x$arguments$exo.var) l <- l-1

    if (max(p, q) > 1) {
        for (i in 2:max(p, q)) {
            H <- cbind(H, H.list[[i]])
        }
        I.n <- diag(ncol(H) - ncol(H.list[[max(p, q)]]))
        Zeros <- matrix(0, nrow = ncol(H) - ncol(H.list[[max(p,
            q)]]), ncol = ncol(H.list[[max(p, q)]]))
        H <- rbind(H, cbind(I.n, Zeros))
        G <- rbind(cbind(G, t(Zeros)), cbind(Zeros, I.n))
        U <- rbind(U, matrix(0, nrow(Zeros), ncol(U)))
    }
    G_inv=solve(G)
    Fmat <- G_inv %*% H
    UtU= cov(t(U))
#    P <- t(chol(t(UtU)))

    # define shock / either in terms of sd or %
    s.j <- rep(0, dim(U)[1])
    ind.j <- rep(0, dim(U)[1])
    
    sigma.il <- vector()
      
    for (i in 1:length(shock.dir))
    {    
      if (shock.var[[i]][1] == 1) {
          j <- shock.var[[i]][2]
      }
      else {
          j <- sum(endoN[1:(shock.var[[i]][1] - 1)]) + shock.var[[i]][2]
      }

      # define variables once
      sigma.il[i] <- sqrt(UtU[j, j])

      if(scal){
#        shock=as.numeric(unlist(strsplit(shock.dir,"/sd"))[[1]])*cons*(1/G_inv[j,j])
        shock <- shock.dir[[i]]*sigma.il[i]/(G_inv%*%UtU)[j,j]
      }
      else{
        shock <- shock.dir[[i]]
      }
      if(!is.numeric(shock)){
        stop("For the argument shock.dir please submit either a volume shock (e.g. shock.dir='0.003/sd') or a standard deviation shock (e.g. shock.dir='-1').")
      }

      s.j[j] <- shock
      ind.j[j] <- 1
    }

    F.n <- list()
    F.n[[1]] <- diag(nrow(Fmat))
    for (i in 2:(n+1)) {
        F.n[[i]] <- F.n[[i-1]] %*% Fmat
    }

    fevd.gi <- matrix(0,length(s.j),n+1)
    fevd.oi <- matrix(0,length(s.j),n+1)
#    for (i in which(s.j==0)[1:(sum(endoN)-1)])
    for (i in 1:sum(endoN))
    {
      e.i <- rep(0,length(s.j))
      e.i[i] <- 1
      temp <- vector()
      temp2 <- vector()
      temp3 <- vector()
      for (j in 1:(n+1))
      {
        temp[j] <- (t(e.i)%*%F.n[[j]]%*%G_inv%*%UtU%*%s.j)^2
        temp2[j] <- t(e.i)%*%F.n[[j]]%*%G_inv%*%UtU%*%t(G_inv)%*%t(F.n[[j]])%*%e.i
#        temp3[j] <- (t(e.i)%*%F.n[[j]]%*%P%*%s.j)^2

        fevd.gi[i,j] <- 1/(sigma.il)*sum(temp[1:j])/sum(temp2[1:j])
#        fevd.oi[i,j] <- 1/(t(e.i)%*%UtU%*%e.i)*sum(temp3[1:j])/sum(temp2[1:j])
      }
    }
    fevd.gi <- fevd.gi[1:dim(x$X)[1],]
    rownames(fevd.gi) <- rownames(x$X)
    fevd.oi <- NULL

    # liste der länder inklusive ländersumme (skaliert, *100 und eine kommastelle)    
    
    gfevd.mat <- t(t(fevd.gi)/apply(fevd.gi,2,sum))
    
    gfevd.list <- list()
    gfevd.list[[1]] <- round(rbind(gfevd.mat[1:endoN[1],],apply(gfevd.mat[1:endoN[1],],2,sum))*100,1)
    colnames(gfevd.list[[1]]) <- 0:n
    rownames(gfevd.list[[1]]) <- c(colnames(x$Data[[1]]),"Sum")
    for (i in 2:l)
    {
    		gfevd.list[[i]] <- round(rbind(gfevd.mat[(sum(endoN[1:(i-1)])+1):(sum(endoN[1:i])),],apply(gfevd.mat[(sum(endoN[1:(i-1)])+1):(sum(endoN[1:i])),],2,sum))*100,1)
    		colnames(gfevd.list[[i]]) <- 0:n
    		rownames(gfevd.list[[i]]) <- c(colnames(x$Data[[i]]),"Sum")
    	}
    names(gfevd.list) <- x$subsys[1:l]
        
    res <- list(gfevd.list = gfevd.list, fevd.gi = fevd.gi, fevd.oi = fevd.oi, Fmat = Fmat, G = G, H = H, U = U,shock.dir=shock.dir,sigma.il=sigma.il)
    return(res)
}

#testFEVD <- FEVD(model1,n=40,list(c(2,1)),shock.dir=list(-1))