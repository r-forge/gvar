# Persitence Profiles (PP) as in Dees, Holly, Pesaran, Smith 2007
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
# pp        ...  matrix of impulse responses
# Fmat      ...  coefficient matrix of ar1 representation of the system, important to check stability of system via its eigenvalues
# G         ...  G
# H         ...  H
# U         ...  U (residual matrix of the gvar system)
# shock.dir ...  contains the submitted shock (either in standard deviations or in other magnitudes that get standardized within GIR)
# sigma.il  ...  residual standard deviation (i.e. magnitude of the shock if not elsewise specified)


PP <- function (x, n=40, shock.var=NULL, shock.dir=NULL, scal=FALSE, system.shock=TRUE){

    endoN=sapply(x$we.vecms,function(x) x$n)
    G <- x$G
    U <- x$U
    W <- x$W
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
        W <- lapply(W,function(y) cbind(y,matrix(0,nrow=nrow(y),ncol=ncol(y))))
    }
    G_inv=solve(G)
    Fmat <- G_inv %*% H
    UtU= cov(t(U))
    
    # define shock / either in terms of sd or %

  if (system.shock!=TRUE)
  {
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
  }

# system wide shock
  if (system.shock==TRUE)
  {
    pp.lm <- vector("list",length=l)
    F.n <- diag(nrow(Fmat))
    covep <- G_inv %*% UtU
    
    for (i in 1:l)
    {
      pp.lm[[i]] <- matrix(0,nrow=x$we.vecm[[i]]$r,ncol=n+1)
      rownames(pp.lm[[i]]) <- paste(1:x$we.vecm[[i]]$r,".",sep="")
      colnames(pp.lm[[i]]) <- 0:n
      for (j in 1:x$we.vecm[[i]]$r)
      {
        if (x$we.vecm[[i]]$r > 1)
        {
          bet <- x$we.vecm[[i]]$beta[,j]
          if (x$arguments$Case[i]=="IV") bet <- x$we.vecm[[i]]$beta[,j][which(names(x$we.vecm[[i]]$beta[,j])!="t")]
          if (x$arguments$Case[i]=="II") bet <- x$we.vecm[[i]]$beta[,j][which(names(x$we.vecm[[i]]$beta[,j])!="constant")]
        } else {
          bet <- as.vector(x$we.vecm[[i]]$beta)
          if (x$arguments$Case[i]=="IV") bet <- x$we.vecm[[i]]$beta[-length(x$we.vecm[[i]]$beta)]
          if (x$arguments$Case[i]=="II") bet <- x$we.vecm[[i]]$beta[-length(x$we.vecm[[i]]$beta)]
        }
        F.n <- diag(nrow(Fmat))
        cons <- bet %*% W[[i]] %*% covep %*% t(W[[i]]) %*% bet
        pp.lm[[i]][j,1] <- (bet %*% W[[i]] %*% F.n %*% covep %*% t(F.n) %*% t(W[[i]]) %*% bet) / cons
        for (step in 1:n)
        {
          F.n <- F.n %*% Fmat
          pp.lm[[i]][j,step+1] <- (bet %*% W[[i]] %*% F.n %*% covep %*% t(F.n) %*% t(W[[i]]) %*% bet) / cons
        } 
      }
    }
    names(pp.lm) <- x$subsys[1:l]
  }

# variable specific shock 

  if (system.shock!=TRUE)
  {
    pp.lm <- vector("list",length=l)
    F.n <- diag(nrow(Fmat))
    covep <- G_inv %*% UtU
    cons <- sigma.il
    
    for (i in 1:l)
    {
      pp.lm[[i]] <- matrix(0,nrow=x$we.vecm[[i]]$r,ncol=n+1)
      rownames(pp.lm[[i]]) <- paste(1:x$we.vecm[[i]]$r,".",sep="")
      colnames(pp.lm[[i]]) <- 0:n  
      for (j in 1:x$we.vecm[[i]]$r)
      {
        if (x$we.vecm[[i]]$r > 1)
        {
          bet <- x$we.vecm[[i]]$beta[,j]
          if (x$arguments$Case[i]=="IV") bet <- x$we.vecm[[i]]$beta[,j][which(names(x$we.vecm[[i]]$beta[,j])!="t")]
          if (x$arguments$Case[i]=="II") bet <- x$we.vecm[[i]]$beta[,j][which(names(x$we.vecm[[i]]$beta[,j])!="constant")]
        } else {
          bet <- as.vector(x$we.vecm[[i]]$beta)
          if (x$arguments$Case[i]=="IV") bet <- x$we.vecm[[i]]$beta[-length(x$we.vecm[[i]]$beta)]
          if (x$arguments$Case[i]=="II") bet <- x$we.vecm[[i]]$beta[-length(x$we.vecm[[i]]$beta)]
        }
        F.n <- diag(nrow(Fmat))
        pp.lm[[i]][j,1] <- (bet %*% W[[i]] %*% F.n %*% covep %*% s.j) / cons
        for (step in 1:n)
        {
          F.n <- F.n %*% Fmat
          pp.lm[[i]][j,step+1] <- (bet %*% W[[i]] %*% F.n %*% covep %*% t(F.n) %*% t(W[[i]]) %*% bet) / cons
        } 
      }
    }
    names(pp.lm) <- x$subsys[1:l]
  }
  
    res <- list(pp = pp.lm, Fmat = Fmat, G = G, H = H, U = U,shock.dir=shock.dir,sigma.il=cons^2)
    return(res)
}



#testPP <- PP(varM_epepsus,n=40)
#testPP <- PP(varM_epepsus,n=40, shock.var=list(c(1,2)), shock.dir=list("1"), scal=FALSE, system.shock=FALSE)
