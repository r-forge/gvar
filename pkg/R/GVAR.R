GVAR <- function(data,tw=NULL,p,q=p,r=NULL,weight,case,exo.var=FALSE,d=NULL,lex=NULL,endo=NULL,ord=NULL,we=NULL,method="max.eigen")
   #     data ... timeseries data as list (each entry is a matrix of a subsystem of variables,
   #              if exo.var=TRUE the last entry are exogeneous variables)
   #       tw ... time window, vector of start and end point
   #        p ... scalar/vector of endogenous lags, (N+1)x1
   #        q ... scalar/vector of weakly exogeneous lags, (N+1)x1
   #        r ... vector of cointegrating relations
   #   weight ... weight matrix of dimension (N+1)x(N+1)
   #     case ... scalar/vector of cases ("I" to "V"), (N+1)x1
   #     endo ... list of endogenous variables used
   #      ord ... vector showing the same variables for weakly exogeneous analysis
   #       we ... list with numbers of weakly exogeneous variables included in each VECM,
   #              corresponds to numbers in ord
   #  exo.var ... if TRUE strictly exogeneous variables are included in the model
   #        d ... list showing which strictly exogeneous variables enter the subsystem equations
   #      lex ... scalar/vector of lags of exogenous variables
   #   method ... select cointegrating rank by max. eigenvalue ("max.eigen") or trace statistic ("trace")

{
  #! ----- Set subsystems -----

  cmodel <- list()

  N <- length(data)-1                                       # number of subsystems i=0,1,...,N

  dims <- vector()
  for (i in 1:(N+1)){dims[i] <- dim(data[[i]])[1]}
  max.dim <- max(dims)
  
  tsi <- tsp(data[[((1:length(dims))[dims==max(dims)])[1]]])
  if (is.null(tw))
  {
    start.ts <- tsi[1]
    end.ts <- tsi[2]
  } else {
    start.ts <- tw[1]
    end.ts <- tw[2]
  }
  freq <- tsi[3]
  

  dt <- 1/freq

  n.exo <- 0
  ex <- 0
  n.ex <- rep(0,N+1)

  if (exo.var)                                              # strictly exogeneous variables
  {
    N <- N-1
    ex <- 1
    n.exo <- dim(data[[length(data)]])[2]
    if (is.null(d))
    {
      d <- vector("list",N+1)
      for (i in 1:(N+1))
      {
        d[[i]] <- 1:n.exo
        n.ex[i] <- n.exo
      }
    } else {
      for (i in 1:(N+1))
      {
      	n.ex[i] <- length(d[[i]])
      }
    }
  }

  if ((length(p)==1)&&(N>=1)){p <- rep(p,(N+1))}
  if ((length(q)==1)&&(N>=1)){q <- rep(q,(N+1))}
  if ((length(case)==1)&&(N>=1)){case <- rep(case,(N+1))}
  if (length(lex)==0){lex <- rep(0,(N+1))} else if ((length(lex)==1)&&(N>=1)){lex <- rep(lex,(N+1))}

  if (is.null(endo))
  {
  	endo <- list()
  	for (i in 1:(N+1))
  	{
  	  endo[[i]] <- 1:dim(data[[i]])[2]
  	}
  }

  n <- vector()                                             # endogeneous variables per subsystem
  for (i in 1:(N+1))
  {
    n[i] <- length(endo[[i]])
  }

  if (is.null(ord))                                         # which variables belong together
  {                                                         # standard: varibales are counted from the beginning of
    ord <- vector()                                         #           each subsystem, last ones are omitted
    for (i in 1:(N+1))
    {
      ord <- c(ord,1:dim(data[[i]])[2])
    }
  }

  if (is.null(we))                                          # which variables enter the VECMs as weakly exogeneous
  {                                                         # standard: all variables that appear in every subsystem
    we <- list()
    fr <- table(ord)
    for (i in 1:(N+1))
    {
      if (max(fr)==(N+1))
      {
        we[[i]] <- as.numeric(names(fr[fr==max(fr)]))
      }else{
        we[[i]] <- vector()
      }
    }
  }
  
  m <- vector()                                             # total variables per subsystem
  for (i in 1:(N+1))
  {
    m[i] <- n[i]+length(we[[i]])
  }

  X <- matrix(nrow=0,ncol=max.dim)                          # create datamatrix X
  for (i in 1:(length(data)-ex))
  {
    X <- rbind(X,t(data[[i]][,endo[[i]]]))
  }

  colnames(X) <- (1-max(p,q)):(max.dim-max(p,q))            # name datamatrix X
  
  if (!is.null(rownames(X)))
  {
    if (is.null(names(data)))
    {
      names(data) <- paste("R",0:N,sep="")
    }
    nam <- vector()
    for (i in 1:(N+1))
    {
      nam <- c(nam,paste(names(data)[i],rownames(X)[(sum(n[1:i])-n[i]+1):sum(n[1:i])],sep="."))
    }
    rownames(X) <- nam
  }else{
    nam <- vector()
    for (i in 0:N)
    {
      for (j in 1:dim(data[[i+1]])[2])
      {
        nam <- c(nam,paste("x",i,".",j,sep=""))
      }
    }
    rownames(X) <- nam
  }
  
  W <- list()                                               # generate weight matrices W_i
  for (i in 1:(N+1))
  {
    W[[i]] <- matrix(0,nrow=m[i],ncol=dim(X)[1])
  }
  for (i in 1:(N+1))
  {
    for (j in 1:n[i])
    {
      if (i==1){l <- 0}else{l <- sum(n[1:i-1])}
      W[[i]][j,l+j] <- 1
    }
    for (j in (n[i]+1):m[i])
    {
      tf <- vector()
      for (l in 1:(N+1))
      {
        tf <- c(tf,is.element(we[[i]][j-n[i]],endo[[l]]))
      }
      W[[i]][j,][ord==we[[i]][j-n[i]]] <- as.numeric(weight[i,tf])
    }
    #W[[i]] <- W[[i]]/apply(W[[i]],1,sum)
  }
  
  #! ----- Calculate VECM models -----
  
  models <- list()
  if (is.null(r)) {r <- vector()}
  rr <- vector()
  for (i in 1:(N+1))
  {
    if (exo.var)
    {
      z <- ts(cbind(t(W[[i]]%*%X),data[[N+2]][,d[[i]]]),start=start.ts,freq=freq)
    } else {
      z <- ts(t(W[[i]]%*%X),start=start.ts,freq=freq)
    }
    cols <- c(colnames(data[[i]])[endo[[i]]],paste(colnames(data[[(1:(N+1))[n==max(n)][1]]]),"*",sep="")[we[[i]]])
    if (exo.var)
    {
      if (!is.null(d[[i]])) {cols <- c(cols,paste(colnames(data[[length(data)]])[d[[i]]],"**",sep=""))}
    }
    if ((n[i]+length(we[[i]])+length(d[[i]]))==length(cols)) {colnames(z) <- cols}
    etw <- list(start=start.ts+max(p,q)*dt,end=end.ts,freq=freq)
    ranks <- rank.test.we(z.ts=z,etw=etw,p=p[i],q=q[i],n=n[i],case=case[i],ex=n.ex[i],lex=lex[i])
    if (is.na(r[i]))
    {
      if (method=="trace")
      {
        for (r_ in 0:(n[i]-1))
        {
          CV <- ranks$CV.trace[[paste("rank",r_,"vs.",n[i])]][1,1]        # niveau set to alpha=5%, for alpha=10%: [1,2], for alpha=5%: [1,3], for alpha=1%: [1,4]
          if ((ranks$LR.trace[[paste("rank",r_,"vs.",n[i])]][1]<CV)&&is.na(r[i])){r[i] <- r_}
        }
      } else {
        for (r_ in 0:(n[i]-1))
        {
          CV <- ranks$CV.maxeigen[[paste("rank",r_,"vs.",r_+1)]][1,1]
          if ((is.na(r[i])&&(ranks$LR.maxeigen[[paste("rank",r_,"vs.",r_+1)]][1]<CV))){r[i] <- r_}
        }
      }
      if (is.na(r[i])){r[i] <- n[i]}
    }
    rr[i] <- r[i]
    if (r[i]==0) {rr[i] <- 1}
    exo <- NULL
    if (!is.null(d[[i]])) {exo=length(d[[i]])}
    mdls <- est.we.mdls(z.ts=z,etw=etw,p=p[i],q=q[i],r=rr[i],n=n[i],case=case[i],ex=n.ex[i],lex=lex[i])
    models[[i]] <- mdls
    names(models)[i] <- names(data)[i]
    cmodel[[i]] <- set.mdl(mdls,exo=exo)
  }
  
  #! ----- Stack to GVAR -----
  
  G <- NULL
  for (i in 1:(N+1))
  {
    Ac <- cbind(diag(n[i]),-cmodel[[i]]$VAR$B_0)
    Wc <- W[[i]]
    G <- rbind(G,Ac%*%Wc)
  }

  H <- vector("list",max(p,q))                                         # different lags allowed, fill with 0's
  for (i in 1:max(p,q))
  {
    for (j in 1:(N+1))
    {
      if (p[j]<i)
      {
        AA <- matrix(0,dim(cmodel[[j]]$VAR$A[[1]])[1],dim(cmodel[[j]]$VAR$A[[1]])[2])
      } else {
        AA <- cmodel[[j]]$VAR$A[[i]]
      }
      if (q[j]<i)
      {
        BB <- matrix(0,dim(cmodel[[j]]$VAR$B[[1]])[1],dim(cmodel[[j]]$VAR$B[[1]])[2])
      } else {
        BB <- cmodel[[j]]$VAR$B[[i]]
      }
      Ac <- cbind(AA,BB)
      H[[i]] <- rbind(H[[i]],Ac%*%W[[j]])
    }
  }

  c.0 <- NULL
  c.1 <- NULL

  if (exo.var)
  {
    Upsilon.0 <- NULL
    Upsilon <- vector("list",max(lex))
  }
  for (i in 1:(N+1))
  {
    if (!is.null(cmodel[[i]]$VAR$c.0))
    {
      c.0 <- rbind(c.0,as.matrix(cmodel[[i]]$VAR$c.0))
    } else c.0<- rbind(c.0,matrix(0,ncol=1,nrow=cmodel[[i]]$VAR$n))
    if (!is.null(cmodel[[i]]$VAR$c.1))
    {
      c.1 <- rbind(c.1,as.matrix(cmodel[[i]]$VAR$c.1))
    } else c.1 <- rbind(c.1,matrix(0,ncol=1,nrow=cmodel[[i]]$VAR$n))
    if (exo.var)                                                           # add coefficient for exogenous variables
    {
      Y.0 <- matrix(0,nrow=n[i],ncol=dim(data[[N+2]])[2])
      if (!is.null(d[[i]])) {Y.0[,d[[i]]] <- as.matrix(cmodel[[i]]$VAR$Upsilon_0)}
      Upsilon.0 <- rbind(Upsilon.0,Y.0)
      for (j in 1:max(lex))
      {
        Y.x <- matrix(0,nrow=n[i],ncol=dim(data[[N+2]])[2])
        if (!is.null(d[[i]]))
        {
          if (lex[i]>=j){Y.x[,1:d[[i]]] <- as.matrix(cmodel[[i]]$VAR$Upsilon[[j]])[,1:d[[i]]]}
        }
        Upsilon[[j]] <- rbind(Upsilon[[j]],Y.x)
      }
    }
  }
  
  T <- dim(X)[2]
  U <- G%*%X[,(max(p,q)+1):T]-matrix(c.0,nrow=dim(X)[1],ncol=T-max(p,q))-c.1%*%t(1:(T-max(p,q)))
  for (k in 1:max(p,q))
  {
    U <- U-H[[k]]%*%X[,(max(p,q)+1-k):(T-k)]
  }
  if (exo.var)
  {
    U <- U-Upsilon.0%*%t(data[[length(data)]][(max(p,q)+1):T,])
    for (k in 1:max(q))
    {
      U <- U-Upsilon[[k]]%*%t(data[[length(data)]][(max(p,q)+1-k):(T-k),])
    }
  }

  U.mean <- apply(U,1,mean)
  U.cov <- tcrossprod(U)/(T-max(p,q))
  
  w.m <- weight/apply(weight,1,sum)
  rownames(w.m) <- colnames(w.m) <- names(data)[1:(N+1)]
  
  if (exo.var==FALSE) {Upsilon.0 <- Upsilon <- NULL}
  
  res <- list(subsys=names(data),G=G,H=H,Upsilon.0=Upsilon.0,Upsilon=Upsilon,c.0=c.0,c.1=c.1,we.vecms=models,m=m,n=n,p=p,q=q,r=rr,weight=w.m,U=U,U.cov=U.cov)
  class(res) <- "GVAR"
  return(res)
}
  