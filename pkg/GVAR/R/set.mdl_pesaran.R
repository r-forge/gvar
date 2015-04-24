set.mdl2 <-function(mdls,exo=NULL,skip=NULL)
# inputs:
# _mdls from a reduced rank VECM estimation
# _chosen rank r
# _number of (strictly) exogenous variables
# outputs:
# _VECM and VAR parameters

{
mdl<- list()
mdl$VECM<- list()
mdl$VAR<- list()
## extract the VECM parameters for the chosen rank r of the cointegration matrix:
for (p in names(mdls)) {
 if (!(p %in% skip)) {
  if (is.list(mdls[[p]])) {
   mdl$VECM[[p]]<- mdls[[p]]
  } else  mdl$VECM[[p]]<- mdls[[p]]
 }
}
mdl$VECM$exo <- exo
                             
## convert the VECM paramters to VAR parameters:
if (mdls[["type"]]=="weakly exogenous VECM") {
 # y_t has n components, x_t has k components, k+n=m
 # VECM: \Delta y_t = c_0+c_1 t+\Lambda\Delta x_t+\sum_{i=1}^{p-1}\Psi_i\Delta z_{t-i}+\Pi.y z_{t-1}+u_t, where u_t is N(0,\Omega_{uu})
 # VAR model:   y_t = c_0+c_1*t+B_0 x_t+\sum_{i=1}^{p}[B_i x_{t-i}+A_i y_{t-i}]+u_t
 m <- mdl$VECM[["m"]]
 n <- mdl$VECM[["n"]]
 ex <- mdl$VECM$ex
 k <- m-n
 p <- mdl$VECM[["p"]]
 q <- mdl$VECM$q
 lex <- mdl$VECM$lex
 
 A <- vector("list",p)
 B <- vector("list",q)
 # NOTE: p>0 is assumed!
 if (p>1) {
  A[[1]]<- diag(n)+mdl$VECM[["Pi.y"]][,1:n]+mdl$VECM[["Phi"]][[1]]
  if (p>2) {
   for (i in 2:(p-1)) {
    A[[i]]<- mdl$VECM[["Phi"]][[i]]-mdl$VECM[["Phi"]][[i-1]]
   }
  }
  A[[p]]<- -mdl$VECM[["Phi"]][[p-1]][,1:n]
 } else if (p==1) {
  A[[1]]<- diag(n)+mdl$VECM[["Pi.y"]][,1:n]
 }
 # NOTE: k>0 is assumed!
 if (q>1) {
  B[[1]]<- -mdl$VECM[["Lambda"]][,1:(m-n)]+mdl$VECM[["Pi.y"]][,(n+1):m]+mdl$VECM[["Psi"]][[1]][,1:k]
  if (q>2) {
   for (i in 2:(q-1)) {
    B[[i]]<- mdl$VECM[["Psi"]][[i]][,1:k]-mdl$VECM[["Psi"]][[i-1]][,1:k]
   }
  }
  B[[q]]<- -mdl$VECM[["Psi"]][[q-1]][,1:k]
 } else if (q==1) {
  B[[1]]<- -mdl$VECM[["Lambda"]][,1:(m-n)]+mdl$VECM[["Pi.y"]][,(n+1):m]
 }
 mdl$VAR$A<- A
 mdl$VAR$B<- B
 mdl$VAR$B_0<- mdl$VECM[["Lambda"]][,1:(m-n)]
 if (!is.null(exo)) {
  Upsilon <- vector("list",lex)
  if ((lex>1) && (q>1)) {
    Upsilon[[1]]<- -mdl$VECM[["Lambda"]][,(m-n+1):(m-n+exo)]+mdl$VECM[["Pi.y"]][,(m+1):(m+exo)]+mdl$VECM[["Psi"]][[1]][,-(1:k)]
    if (lex>2) {
      for (i in 2:(lex-1)) {
        if (q>=i)
        {
          Upsilon[[i]] <- mdl$VECM[["Psi"]][[i]][,-(1:k)]-mdl$VECM[["Psi"]][[i-1]][,-(1:k)]
        } else {
          Upsilon[[i]] <- mdl$VECM[["Psi"]][[i]][,1:exo]-mdl$VECM[["Psi"]][[i-1]][,1:exo]
        }
      }
    }
    if (q>=lex) {
      Upsilon[[lex]] <- -mdl$VECM[["Psi"]][[lex-1]][,-(1:k)]
    } else {
      Upsilon[[lex]] <- -mdl$VECM[["Psi"]][[lex-1]][,1:exo]
    }
  } else if ((lex>1) && (q<2)) {
    Upsilon[[1]] <- -mdl$VECM[["Lambda"]][,1:exo]+mdl$VECM[["Pi.y"]][,(m+1):(m+exo)]+mdl$VECM[["Psi"]][[1]][,1:exo]
    if (lex>2) {
      for (i in 2:(lex-1)) {
        Upsilon[[i]] <- mdl$VECM[["Psi"]][[i]][,1:exo]-mdl$VECM[["Psi"]][[i-1]][,1:exo]
      }      
    }
    Upsilon[[lex]] <- -mdl$VECM[["Psi"]][[lex-1]][,1:exo]
  } else if (lex==1) {
    Upsilon[[1]]<- -mdl$VECM[["Lambda"]][,(m-n+1):(m-n+exo)]+mdl$VECM[["Pi.y"]][,(m+1):(m+exo)]
  }
  mdl$VAR$Upsilon<- Upsilon
  mdl$VAR$Upsilon_0 <- mdl$VECM[["Lambda"]][,(m-n+1):(m-n+exo)]
 }

 for (pn in names(mdl$VECM)) {
  if (!(pn %in% c("Lambda","Psi","alpha","beta","Pi.y"))) {
   mdl$VAR[[pn]]<- mdl$VECM[[pn]]
  }
 }
 mdl$VAR$exo <- exo
} else if ( mdls[["type"]]=="pure VECM" ) {
 # Y_t has n components
 # VECM: \Delta Y_t = \Pi Y_{t-1}+\sum_{i=1}^{k-1} \Gamma_i\Delta Y_{t-i}+\Phi D_t+\epsilon_t, where \epsilon_t is N(0,\Omega) and \Phi D_t= \mu_0+\mu_1 t
 # VAR model:   Y_t = A_1 Y_{t-1}+...+ A_k Y_{t-k} +\Phi D_t+\epsilon_t
 p<- mdl$VECM[["p"]]
 n<- mdl$VECM[["n"]]
 A<- vector("list",p)
 # NOTE: k>0 is assumed!
 if (p>1) {
  A[[1]]<- diag(n)+mdl$VECM[["Pi"]]+mdl$VECM[["Gamma"]][[1]]
  if (p>2) {
   for (i in 2:(p-1)) {
    A[[i]]<- mdl$VECM[["Gamma"]][[i]]-mdl$VECM[["Gamma"]][[i-1]]
   }
  }
  A[[p]]<- -mdl$VECM[["Gamma"]][[p-1]]
 } else if (p==1) {
  A[[1]]<- diag(n)+mdl$VECM[["Pi"]]
 }
 mdl$VAR$A<- A
 for (pn in names(mdl$VECM)) {
  if (!any(pn %in% c("Gamma","alpha","beta","Pi"))) {
   mdl$VAR[[pn]]<- mdl$VECM[[pn]]
  }
 }
}

return(mdl)

}