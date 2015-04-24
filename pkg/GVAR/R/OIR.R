# Orthogonalized Impulse Response Function (GIR) as in Peseran and Shin 1996
#
# input
# x         ...  a gvar object 
# n         ...  number of periods, in the GVAR literature often set to 40  
# shock.var ...  list of vectors of length 2, first element specifies the country to be shocked, second element the variable in the country (e.g. shock.var=c(13,3) third variable in US is going to be shocked)
# shock.dir ...  a list of scalars, set to 1 or -1, for a (negative) 1sd shock, where sd refers to the country specific residuals. Alternatively, if scal = TRUE, any number submitted indicates 
#                the magnitude of the resulting first period shock, e.g. -0.01 for a 1% decrease in the first quarter
# scal      ...  TRUE if shock.dir indicates resulting value of shock
#
# output
# psi       ...  matrix of impulse responses
# Fmat      ...  coefficient matrix of ar1 representation of the system, important to check stability of system via its eigenvalues
# G         ...  G
# H         ...  H
# U         ...  U (residual matrix of the gvar system)
# shock.dir ...  contains the submitted shock (either in standard deviations or in other magnitudes that get standardized within GIR)
# sigma.il  ...  residual standard deviation (i.e. magnitude of the shock if not elsewise specified)

# load("Z:/Documents/GVAR/martin/20110309/varM_epepsus.rda")

OIR=function (x, n=40, shock.var, shock.dir=-1,scal=FALSE){

    shock.var <- list(c(1,shock.var))
    shock.dir <- list(shock.dir)

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
    
    U.0 <- t(U[1:endoN[1],])
    covU.0 <- cov(U.0)
    
    P.0 <- t(chol(t(covU.0)))
    
    P.0G <- diag(dim(H)[1])
    P.0G[1:endoN[1],1:endoN[1]] <- P.0
    
    V <- P.0G%*%U
    covV <- cov(t(V)) 
    
    #P <- t(chol(t(UtU)))
    
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
    
    P.0G_inv <- solve(P.0G)
    cons <- sqrt(as.vector(s.j%*%covV%*%s.j))

    psi.m <- matrix(nrow = nrow(U), ncol = n + 1)
    F.n <- diag(nrow(Fmat))
    psi.m[,1] <- (F.n %*% G_inv %*% P.0G_inv %*% covV %*% s.j)/cons
    for (i in 1:n)
    {
      F.n <- F.n %*% Fmat
      psi.m[,(i+1)] <- (F.n %*% G_inv %*% P.0G_inv %*% covV %*% s.j)/cons
    }
    
    psi <- list()
    psi[[1]] <- t(psi.m[1:endoN[1],])
    rownames(psi[[1]]) <- 0:n
    colnames(psi[[1]]) <- colnames(x$data[[1]])
    for (i in 2:l)
    {
    		psi[[i]] <- t(psi.m[(sum(endoN[1:(i-1)])+1):(sum(endoN[1:i])),])
    		rownames(psi[[i]]) <- 0:n
    		colnames(psi[[i]]) <- colnames(x$Data[[i]])  
    	}
    names(psi) <- x$subsys[1:l]
    
    res <- list(psi = psi, Fmat = Fmat, G = G, H = H, U = U,shock.dir=shock.dir,sigma.il=sigma.il)
    return(res)
}

#testOIR <- OIR(varM_epepsus,n=40,shock.var=1,shock.dir=-1,scal=FALSE) 