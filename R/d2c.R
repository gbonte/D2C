#' @import MASS randomForest corpcor xgboost lazy




npred<-function(X,Y,lin=TRUE,norm=TRUE){
  ## normalized mean squared error of the dependency
  N<-NROW(X)
  n<-NCOL(X)
  
  if (n>1){
    w.const<-which(apply(X,2,sd)<0.01)    
    if (length(w.const)>0){
      X<-X[,-w.const]
    }
    n<-NCOL(X)
  }
  
  if (lin)
    return(max(1e-3,regrlin(X,Y)$MSE.loo/var(Y)))
  
  
  XX<-scale(X)
  if (N<5 | any(is.na(XX)))
    stop("Error in npred")
  e<-Y-lazy.pred(XX,Y,XX,conPar=c(min(10,N-2),min(N,20)),
                 linPar=NULL,class=FALSE,cmbPar=10)
  #  Itr=sample(1:N,min(50,round(2*N/3)))
  #  Its=setdiff(1:N,Itr)
  # e<-Y[Its]-rf.pred(X[Itr,],Y[Itr],X[Its,],class=FALSE,ntree=5)
  if (norm){
    nmse<-mean(e^2)/var(Y) 
    return(max(1e-3,nmse))
  }
  
  return(mean(e^2))
}


norminf<-function(y,x1,x2=NULL,lin=TRUE){
  ## Normalized conditional information of x1 to y given x2
  ## I(x1;y| x2)= (H(y|x2)-H(y | x1,x2))/H(y|x2)
  
  if (is.null(x2)){ ## I(x1;y)= (H(y)-H(y | x1))/H(y)
    return(max(0,1-npred(x1,y,lin=lin,norm=TRUE)))
    
  }
  
  np<-npred(x2,y,lin=lin,norm=FALSE)
  x1x2<-cbind(x1,x2)
  delta<- max(0,np-npred(x1x2,y,lin=lin,norm=FALSE))/(np+0.01)
  return(delta)
  
}


#' compute descriptor
#' @param D :  the observed data matrix of size [N,n], where N is the number of samples and n is the number of nodes
#' @param ca : node index (\eqn{1 \le ca \le n}) of the putative cause
#' @param ef : node index (\eqn{1 \le ef \le n}) of the putative effect
#' @param ns : size of the Markov Blanket
#' @param lin : TRUE OR FALSE. if TRUE it uses a linear model to assess a dependency, otherwise a local learning algorithm 
#' @param acc : TRUE OR FALSE. if TRUE it uses the accuracy of the regression as a descriptor
#' @param struct :   TRUE or FALSE to use the ranking in the markov blanket as a descriptor
#' @param pq :  a vector of quantiles used to compute de descriptor
#' @param bivariate :  TRUE OR FALSE. if TRUE it includes the descriptors of the bivariate dependency
#' @param maxs : max number of pairs MB(i), MB(j) considered 
#' @param boot :  feature selection algorithm
#' @details This function is the core of the D2C algorithm. Given two candidate nodes, (\code{ca}, putative cause and \code{ef}, putative effect) it first infers from the dataset D the Markov Blankets of the variables indexed by \code{ca} and \code{ef} (\code{MBca} and \code{MBef}) by using the \link{mimr} algorithm (Bontempi, Meyer, ICML10). Then it computes a set of (conditional) mutual information terms describing the dependency between the variables ca and ef. These terms are used to create a vector of descriptors. If \code{acc=TRUE}, the vector contains the descriptors related to the asymmetric information theoretic terms described in the paper. If \code{struct=TRUE}, the vector contains descriptors related to the positions of the terms of the MBef in MBca and viceversa. The estimation of the information theoretic terms require the estimation of the dependency between nodes. If \code{lin=TRUE} a linear assumption is made. Otherwise the local learning estimator, implemented by the R package \link{lazy}, is used.
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @references Bontempi G., Meyer P.E. (2010) Causal filter selection in microarray data. ICML'10
#' @references M. Birattari, G. Bontempi, and H. Bersini (1999) Lazy learning meets the recursive least squares algorithm. Advances in Neural Information Processing Systems 11, pp. 375-381. MIT Press.
#' @references G. Bontempi, M. Birattari, and H. Bersini (1999) Lazy learning for modeling and control design. International Journal of Control, 72(7/8), pp. 643-658.
#' @export 
descriptor<-function(D,ca,ef,ns=min(4,NCOL(D)-2),
                     lin=FALSE,acc=TRUE,struct=FALSE, 
                     pq= c(0.1,0.25,0.5,0.75,0.9),
                     bivariate=FALSE,maxs=10,boot="mimr" ){
  if (bivariate){
    return(c(D2C.n(D,ca,ef,ns,lin,acc,struct,pq=pq,boot=boot,maxs=maxs),D2C.2(D[,ca],D[,ef])))
  }else {
    De=D2C.n(D,ca,ef,ns,lin,acc,struct,pq=pq,boot=boot,maxs=maxs)
    De=c(NROW(D),NCOL(D)/NROW(D),kurtosis(D[,ca])/kurtosis(D[,ef]),kurtosis(D[,ef])/kurtosis(D[,ca]),De)
    if (length(names(De))<=4){
      print(De)
      stop("Error in descriptor")
    }
    names(De)[1:4]=c('N', 'n/N','kurtosis1','kurtosis2')
    return(De)
    }
}




E<-function(x){
  return(ecdf(x)(x))
}

D2C.n<-function(D,ca,ef,ns=min(4,NCOL(D)-2),maxs=20,
                lin=FALSE,acc=TRUE,struct=TRUE,
                pq= c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),boot="mrmr"){
  ## is i cause oj j
  n<-NCOL(D)
  N<-NROW(D)
  
  #### creation of the Markov Blanket of ca (denoted MBca)
  #### MB is obtained by first ranking the other nodes and then selecting a subset of size ns 
  #### with the algorithm mentioned in the "boot" variable
  
  
  MBca<-setdiff(1:n,ca)
  MBef<-setdiff(1:n,ef)
  MBca2=MBca
  MBef2<-MBef
  if (n>(ns+1)){
    ind<-setdiff(1:n,ca)
    ind<-ind[rankrho(D[,ind],D[,ca],nmax=min(length(ind),5*ns))]
    if (boot=="mrmr")
      MBca<-ind[mrmr(D[,ind],D[,ca],nmax=ns)]
    
    if (boot=="rank")
      MBca<-ind[1:ns]
    
    MBca2=MBca
    
    if (boot=="mimr"){
      MBca2<-ind[mimreff(D[,ind],D[,ca],nmax=ns)] ## putative list of effects
      MBca<-ind[mimr(D[,ind],D[,ca],nmax=2*ns,init=TRUE)] 
      MBca<-setdiff(MBca,MBca2) ## remove putative effects
      MBca<-MBca[1:min(length(MBca),ns)]
    }
    
    
    
    #### creation of the Markov Blanket of ef (denoted MBef)
    ind2<-setdiff(1:n,ef)
    ind2<-ind2[rankrho(D[,ind2],D[,ef],nmax=min(length(ind2),5*ns))]
    
    if (boot=="mrmr")
      MBef<-ind2[mrmr(D[,ind2],D[,ef],nmax=ns)]
    
    if (boot=="rank")
      MBef<-ind2[1:ns]
    
    MBef2<-MBef
    if (boot=="mimr"){
      MBef2<-ind2[mimreff(D[,ind2],D[,ef],nmax=ns)]
      MBef<-ind2[mimr(D[,ind2],D[,ef],nmax=2*ns,init=TRUE)]
      MBef<-setdiff(MBef,MBef2)   ## remove putative effects
      MBef<-MBef[1:min(ns,length(MBef))]
    }
    
    if (boot=="mrmr2"){
      ind<-setdiff(union(ind,ind2),c(ca,ef))
      fs=mrmr2(D[,ind],D[,ca],D[,ef],nmax=ns)
      MBca<-ind[fs$fs1] ## putative list of effects
      MBef<-ind[fs$fs2] 
    }
    
    
    iMb=intersect(MBca,MBef)

    if (length(iMb)<(length(MBca)-2))
      MBca=setdiff(MBca,iMb)
    
    if (length(iMb)<(length(MBef)-2))
      MBef=setdiff(MBef,iMb)
    
   if (length(MBca)<3 || length(MBef)<3)
     stop("error: too short MB size")
  }
  
  namesx<-NULL 
  x<-NULL
  
  
  
  if (struct){
    
    ## position of effect in the MBca 
    if (is.element(ef,MBca))
      pos.ef<-(which(MBca==ef))/(ns)
    else
      pos.ef<-2
    
    ## position of ca in the MBef
    if (is.element(ca,MBef))
      pos.ca<-(which(MBef==ca))/(ns)
    else
      pos.ca<-2
    
    sx.ef<-NULL
    ## position of variables of MBef in MBca 
    for (i in 1:length(MBef)){
      if (is.element(MBef[i],MBca))
        sx.ef<-c(sx.ef,(which(MBca==MBef[i]))/(ns))
      else
        sx.ef<-c(sx.ef,2)
      
    }
    
    ## position of variables of MBca in MBef
    sx.ca<-NULL
    for (i in 1:length(MBca)){
      if (is.element(MBca[i],MBef))
        sx.ca<-c(sx.ca,(which(MBef==MBca[i]))/(ns))
      else
        sx.ca<-c(sx.ca,2)      
    }
    
    
    x<-c(x,pos.ca,pos.ef,quantile(sx.ca,probs=pq),quantile(sx.ef,probs=pq))
    namesx<-c(namesx,"pos.ca","pos.ef",paste0("sx.ca",1:length(pq)),
              paste0("sx.ef",1:length(pq)))
  } ## if struct
  
  MBca<-setdiff(MBca,ef)
  MBef<-setdiff(MBef,ca)
  
  MBca2<-setdiff(MBca2,ef)
  MBef2<-setdiff(MBef2,ca)
  
  if (acc){
    
    ## relevance of ca for ef
    ca.ef<-norminf(D[,ca],D[,ef],lin=lin) #I(zi;zj)
    ## relevance of ef for ca
    ef.ca<-norminf(D[,ef],D[,ca],lin=lin) #I(zj;zi)
    
    
    delta<- norminf(D[,ef],D[,ca],D[,MBef],lin=lin) #I(zi;zj|Mj)
    delta2<- norminf(D[,ca],D[,ef],D[,MBca],lin=lin) #I(zi;zj|Mi)
    
    #Delta<- norminf(D[,ef],D[,ca],D[,MBca],lin=lin) #I(zi;zj|Mj)
    #Delta2<- norminf(D[,ca],D[,ef],D[,MBef],lin=lin) #I(zi;zj|Mi)
    
    ## relevance of ca for ef given MBef
    delta.i<-NULL
    for (m in MBef)
      delta.i<- c(delta.i,norminf(D[,ef],D[,ca],D[,c(m,MBca)],lin=lin)) #I(zj;zi|Mj^k)
    
    delta2.i<-NULL
    for (m in MBca)
      delta2.i<- c(delta2.i,norminf(D[,ca],D[,ef],D[,c(m,MBef)],lin=lin)) #I(zi;zj|Mi^k)
    
    
    
    I1.i<-NULL
    ## Information of Mbef on ca (i)
    
    for (j in 1:length(MBef)){
      I1.i<-c(I1.i, (norminf(D[,MBef[j]],D[,ca],lin=lin)))  ## I(Mj^k;zi) equation (11)
    }
    
    I1.j<-NULL
    ## Information of Mbca on ef  (j)
    
    for (j in 1:length(MBca)){
      I1.j<-c(I1.j, (norminf(D[,MBca[j]],D[,ef],lin=lin))) ## I(Mi^k;zj) equation (11)
    }
    
    I2.i<-NULL
    I2.ib<-NULL
    ## Information of Mbef on ca given ef
    for (j in 1:length(MBef)){
      backdoor=setdiff(union(MBca,MBef[-j]),c(ca,ef))
      I2.i<-c(I2.i, norminf(D[,ca], D[,MBef[j]],D[,ef],lin=lin)) ## I(zi; Mj^k|zj) equation (8)
      I2.ib<-c(I2.ib, norminf(D[,ca], D[,MBef[j]],D[,c(ef,backdoor)],lin=lin)) 
    }
    
    I2.j<-NULL
    I2.jb<-NULL
    ## Information of Mbca on ef given ca
    for (j in 1:length(MBca)){
      backdoor=setdiff(union(MBef,MBca[-j]),c(ca,ef))
      I2.j<-c(I2.j, norminf(D[,ef], D[,MBca[j]],D[,ca],lin=lin)) ## I(zj; Mi^k|zi) equation (8)
      I2.jb<-c(I2.jb, norminf(D[,ef], D[,MBca[j]],D[,c(ca,backdoor)],lin=lin)) 
    }
    
    IJ<-expand.grid(1:length(MBca),1:length(MBef))
    IJ<-IJ[sample(1:NROW(IJ),min(maxs,NROW(IJ))),]
    
    I3.i<-NULL
    I3.ib<-NULL
    ## Information of MBef on MBca given ca
    for (r in 1:NROW(IJ)){
      i=IJ[r,1]
      j=IJ[r,2]
      backdoor=setdiff(union(MBca[-i],MBef[-j]),ca)
      I3.ib<-c(I3.ib,(norminf(D[,MBca[i]],D[,MBef[j]],D[,c(ca,backdoor)],lin=lin))) 
      I3.i<-c(I3.i,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ca],lin=lin))) ## I(Mi^k; Mj^k|zi) equation (9-10)
      
    }
    
    I3.j<-NULL
    I3.jb<-NULL
    ## Information of MBef on MBca given ef
    for (r in 1:NROW(IJ)){
      i=IJ[r,1]
      j=IJ[r,2]
      backdoor=setdiff(union(MBca[-i],MBef[-j]),ef)
      I3.jb<-c(I3.jb,(norminf(D[,MBca[i]],D[,MBef[j]],D[,c(ef,backdoor)],lin=lin))) 
      I3.j<-c(I3.j,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ef],lin=lin))) ## I(Mi^k; Mj^k|zj) equation (9-10)
      
      
    }
    
    
    IJ<-expand.grid(1:length(MBca),1:length(MBca))
    IJ<-IJ[which(IJ[,1]<IJ[,2]),]
    IJ<-IJ[sample(1:NROW(IJ),min(maxs,NROW(IJ))),]
    
    Int3.i<-NULL
    ## Interaction of terms of Mbca
    for (r in 1:NROW(IJ)){
      i=IJ[r,1]
      j=IJ[r,2]
      Int3.i<-c(Int3.i,(norminf(D[,MBca[i]],D[,MBca[j]],D[,ca],lin=lin)
                        -norminf(D[,MBca[i]],D[,MBca[j]],lin=lin))) ## I(Mi^k; Mi^k|zi)-I(Mi^k; Mi^k)
    }
    
    IJ<-expand.grid(1:length(MBef),1:length(MBef))
    IJ<-IJ[which(IJ[,1]<IJ[,2]),]
    IJ<-IJ[sample(1:NROW(IJ),min(maxs,NROW(IJ))),]
    
    Int3.j<-NULL
    ## Interaction of terms of Mbef
    for (r in 1:NROW(IJ)){
      i=IJ[r,1]
      j=IJ[r,2]
      Int3.j<-c(Int3.j,(norminf(D[,MBef[i]],D[,MBef[j]],D[,ef],lin=lin)
                        -norminf(D[,MBef[i]],D[,MBef[j]],lin=lin))) ## I(Mj^k; Mj^k|zj)-I(Mj^k; Mj^k)
    }
    
    E.ef=ecdf(D[,ef])(D[,ef]) ## empirical cdf of D[,ef]
    E.ca=ecdf(D[,ca])(D[,ca])
    
    ## gini relevance of ca for ef
    gini.ca.ef<-norminf(D[,ca],E.ef,lin=lin)
    ## gini relevance of ef for ca
    gini.ef.ca<-norminf(D[,ef],E.ca,lin=lin)
    gini.delta<- norminf(D[,ef],E.ca,D[,MBef],lin=lin)
    gini.delta2<- norminf(D[,ca],E.ef,D[,MBca],lin=lin)
    
    
    if (FALSE){
      
      
      Int1.i<-NULL
      ## 
      for (j in 1:length(MBef)){
        Int1.i<-c(Int1.i, norminf(D[,ca], D[,MBef[j]],D[,ef],lin=lin)-norminf(D[,ca], D[,MBef[j]],lin=lin)) ## I(zi; Mj^k|zj)-I(zi; Mj^k)
      }
      
      Int1.j<-NULL
      ## 
      for (j in 1:length(MBca)){
        Int1.j<-c(Int1.j, norminf(D[,ef], D[,MBca[j]],D[,ca],lin=lin)
                  -norminf(D[,ef], D[,MBca[j]],lin=lin)) ## I(zj; Mi^k|zi)- I(zj; Mi^k)
      }
      
      
      Int2.i<-NULL
      ## 
      for (r in 1:NROW(IJ)){
        i=IJ[r,1]
        j=IJ[r,2]
        Int2.i<-c(Int2.i,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ca],lin=lin)
                          -norminf(D[,MBca[i]],D[,MBef[j]],lin=lin))) ## I(Mi^k; Mj^k|zi)-I(Mi^k; Mj^k)
      }
      
      Int2.j<-NULL
      ## 
      for (r in 1:NROW(IJ)){
        i=IJ[r,1]
        j=IJ[r,2]
        Int2.j<-c(Int2.j,(norminf(D[,MBca[i]],D[,MBef[j]],D[,ef],lin=lin)
                          -norminf(D[,MBca[i]],D[,MBef[j]],lin=lin))) ## I(Mi^k; Mj^k|zj)-I(Mi^k; Mj^k)
      }
      
      
      E3.i<-NULL
      ## Information of MBef2 on MBca given ca
      for (i in 1:length(MBca))
        for (j in 1:length(MBef2)){
          E3.i<-c(E3.i,(norminf(D[,MBca[i]],D[,MBef2[j]],D[,ca],lin=lin)))
        }
      
      E3.j<-NULL
      ## Information of MBef on MBca2 given ef
      for (i in 1:length(MBca2))
        for (j in 1:length(MBef)){
          E3.j<-c(E3.j,(norminf(D[,MBca2[i]],D[,MBef[j]],D[,ef],lin=lin)))
        }
      
      
      G1.i<-NULL
      ## Information of Mbef on ca 
      
      for (j in 1:length(MBef)){
        G1.i<-c(G1.i, (npred(D[,MBef[j]],E(D[,ca]),lin=lin)))
      }
      
      G1.j<-NULL
      ## Information of Mbca on ef 
      
      for (j in 1:length(MBca)){
        G1.j<-c(G1.j, (npred(D[,MBca[j]],E(D[,ef]),lin=lin)))
      }
      
      G2.i<-NULL
      ## Information of Mbef on ca given ef
      for (j in 1:length(MBef)){
        G2.i<-c(G2.i, norminf(D[,ca], E(D[,MBef[j]]),D[,ef],lin=lin))
      }
      
      G2.j<-NULL
      ## Information of Mbca on ef given ca
      for (j in 1:length(MBca)){
        G2.j<-c(G2.j, norminf(D[,ef], E(D[,MBca[j]]),D[,ca],lin=lin))
      }
      
      
      G3.i<-NULL
      ## Information of MBef on MBca given ca
      for (i in 1:length(MBca))
        for (j in 1:length(MBef)){
          G3.i<-c(G3.i,(norminf(D[,MBca[i]],E(D[,MBef[j]]),D[,ca],lin=lin)))
        }
      
      G3.j<-NULL
      ## Information of MBef on MBca given ef
      for (i in 1:length(MBca))
        for (j in 1:length(MBef)){
          G3.j<-c(G3.j,(norminf(D[,MBca[i]],E(D[,MBef[j]]),D[,ef],lin=lin)))
        }
    } ## if FALSE
    
    x<-c(x,delta,delta2,
         quantile(delta.i,probs=pq,na.rm=TRUE),
         quantile(delta2.i,probs=pq,na.rm=TRUE),ca.ef,ef.ca,
         quantile(I1.i,probs=pq,na.rm=TRUE),quantile(I1.j,probs=pq,na.rm=TRUE),
         quantile(I2.i,probs=pq,na.rm=TRUE),quantile(I2.j,probs=pq,na.rm=TRUE),
         quantile(I2.ib,probs=pq,na.rm=TRUE),quantile(I2.jb,probs=pq,na.rm=TRUE),
         quantile(I3.i,probs=pq,na.rm=TRUE),quantile(I3.j,probs=pq,na.rm=TRUE),
         quantile(I3.ib,probs=pq,na.rm=TRUE),quantile(I3.jb,probs=pq,na.rm=TRUE),
         quantile(Int3.i,probs=pq,na.rm=TRUE),quantile(Int3.j,probs=pq,na.rm=TRUE),
         gini.delta,gini.delta2,
         gini.ca.ef,gini.ef.ca
         #quantile(Int1.i,probs=pq,na.rm=TRUE),quantile(Int1.j,probs=pq,na.rm=TRUE),
         #quantile(Int2.i,probs=pq,na.rm=TRUE),quantile(Int2.j,probs=pq,na.rm=TRUE),
         # quantile(G1.i,probs=pq,na.rm=TRUE),quantile(G1.j,probs=pq,na.rm=TRUE),
         # quantile(G2.i,probs=pq,na.rm=TRUE),quantile(G2.j,probs=pq,na.rm=TRUE),
         # quantile(G3.i,probs=pq,na.rm=TRUE),quantile(G3.j,probs=pq,na.rm=TRUE),
         #quantile(E3.i,probs=pq,na.rm=TRUE),quantile(E3.j,probs=pq,na.rm=TRUE)
    )
    
    namesx<-c(namesx,"delta","delta2",
              paste("delta.i",1:length(pq)),
              paste("delta2.i",1:length(pq)),
              "ca.ef","ef.ca",
              paste0("I1.i",1:length(pq)), paste0("I1.j",1:length(pq)),
              paste0("I2.i",1:length(pq)), paste0("I2.j",1:length(pq)),
              paste0("I2.ib",1:length(pq)), paste0("I2.jb",1:length(pq)),
              paste0("I3.i",1:length(pq)), paste0("I3.j",1:length(pq)),
              paste0("I3.ib",1:length(pq)), paste0("I3.jb",1:length(pq)),
              paste0("Int3.i",1:length(pq)), paste0("Int3.j",1:length(pq)),
              "gini.delta","gini.delta2",
              "gini.ca.ef","gini.ef.ca"
              #paste0("Int1.i",1:length(pq)), paste0("Int1.j",1:length(pq)),
              #paste0("Int2.i",1:length(pq)), paste0("Int2.j",1:length(pq)),
              #   paste0("G1.i",1:length(pq)), paste0("G1.j",1:length(pq)),
              #    paste0("G2.i",1:length(pq)), paste0("G2.j",1:length(pq)),
              #    paste0("G3.i",1:length(pq)), paste0("G3.j",1:length(pq)),
              #paste0("E3.i",1:length(pq)), paste0("E3.j",1:length(pq))
    )
  } ## if acc
  
  
  if (length(namesx)!=length(x)){
    print(x)
    stop("error in D2C.n")
  }
  names(x)<-namesx
  x
  
  
}








D2C.2<-function(x,y,sdkmin=0.5,sdkmax=0.5,Ls=1){
  
  
  
  copula<-TRUE
  alpha<-FALSE
  
  origx<-x
  origy<-y
  
  x<-scale(x)
  y<-scale(y)
  
  sdk<-mean(c(sdkmin,sdkmax))
  
  ##############################@
  N<-length(x)
  point.per.breaks<-min(150,round(N/3))
  if (alpha){
    ncp<-max(2,sqrt(N/point.per.breaks))
    ind.cpx<-seq(min(x),1.01*max(x),length.out=ncp+1)
    ind.cpy<-seq(min(y),1.01*max(y),length.out=ncp+1)
    
    cp<-array(NA,c(ncp,ncp))
    for (i in 1:(ncp))
      for (j in 1:(ncp)){
        
        cp[i,j]<-length(which(x<ind.cpx[i+1] & x>=ind.cpx[i] & y<ind.cpy[j+1] & y>=ind.cpy[j]))/N
        
      }
    
    
    
    alpha<-array(NA,c(ncp-1,ncp-1))
    alpha2<-array(NA,c(ncp-1,ncp-1))
    for (i in 1:(ncp-1))
      for (j in 1:(ncp-1)){
        alpha[i,j]<-log(cp[i,j]*cp[i+1,j+1]/(cp[i,j+1]*cp[i+1,j]))
        
      }
    alpha[which(is.nan(alpha))]<-0
    alpha[which(is.infinite(alpha))]<-0
    alpha[which(alpha==0)]<-NA
    qalpha<-quantile(c(alpha),qnt,na.rm=TRUE)
  }
  
  
  ###################### inp: x1, out: y1, err: Ey
  
  kp<-150
  ix<-sort(x,decreasing=FALSE,ind=TRUE)$ix
  x1<-x[ix]
  y1<-y[ix]
  
  Hy<-NULL
  Ey<-NULL
  ypred<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    
    if (Ls>1){
      S<-seq(sdkmin,sdkmax,length.out=Ls)
      el<-NULL
      for (ssdk in S){
        wd<-numeric(N)
        wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*ssdk^2))
        if (sum(wd)>0.01){
          wd<-wd/sum(wd)
          py<-sum(wd*y1)
          el<-c(el,abs(py-y1[i]))
        } else {
          el<-c(el,Inf)
        }
      }
      
      ssdk<-S[which.min(el)]
    } else {
      ssdk<-sdk
    }
    wd<-numeric(N)
    wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*ssdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      py<-sum(wd*y1)
      sdy<-sqrt(sum(wd*((y1-py)^2)))
    } else {
      py<-mean(y1)
      sdy<-1
    }
    Hy<-c(Hy,sd(y1)-sdy)
    Ey<-c(Ey,y1[i]-py)
    ypred<-c(ypred,py)
  }
  
  EEy<-NULL
  HHy<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    wd<-numeric(N)
    wd[ind]<- exp(-((x1[ind]-x1[i])^2)/(2*sdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      
      py<-sum(wd*Ey)
      sdy<-sqrt(sum(wd*((Ey-py)^2)))
    } else {
      py<-mean(Ey)
      sdy<-1
    }
    HHy<-c(HHy,sd(Ey)-sdy)
    
    EEy<-c(EEy,Ey[i]-py)
    
  }
  
  
  ##################################@ inp y2: out x2, err: Ex 
  iy<-sort(y,decreasing=FALSE,ind=TRUE)$ix
  x2<-x[iy]
  y2<-y[iy]
  
  Ex<-NULL
  Hx<-NULL
  xpred<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    if (Ls>1){
      S<-seq(sdkmin,sdkmax,length.out=Ls)
      el<-NULL
      for (ssdk in S){
        wd<-numeric(N)
        wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*ssdk^2))
        
        if (sum(wd)>0.01){
          wd<-wd/sum(wd)
          px<-sum(wd*x2)
          el<-c(el,abs(px-x2[i]))
        } else {
          el<-c(el,Inf)
        }
      }
      
      ssdk<-S[which.min(el)]
    } else {
      ssdk<-sdk
    }
    
    wd<-numeric(N)
    wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*ssdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      px<-sum(wd*x2)
      sdx<-sqrt(sum(wd*((x2-px)^2)))
    } else {
      px<-x2[i]
      sdx<-1
    }
    
    Hx<-c(Hx,sd(x2)-sdx)    
    Ex<-c(Ex,x2[i]-px)
    xpred<-c(xpred,px)
    
  }
  
  EEx<-NULL
  HHx<-NULL
  for (i in 1:N ){
    ind<-setdiff(max(1,(i-kp)):min(N,(i+kp)),i)
    wd<-numeric(N)
    wd[ind]<- exp(-((y2[ind]-y2[i])^2)/(2*sdk^2))
    if (sum(wd)>0){
      wd<-wd/sum(wd)
      px<-sum(wd*Ex)
    } else {
      px<-Ex[i]
      sdx<-1
      
    }
    sdx<-sqrt(sum(wd*((Ex-px)^2)))
    HHx<-c(HHx,sd(Ex)-sdx)
    EEx<-c(EEx,Ex[i]-px)
    
  }
  
  mx<-c(median(Hx),max(Hx))
  dx<-c(max(Hx)-min(Hx))
  
  my<-c(median(Hy),max(Hy))
  dy<-c(max(Hy)-min(Hy))
  
  qnt<-seq(0.05,0.99,length=3)
  
  qx<-quantile(origx,qnt)
  qy<-quantile(origy,qnt)
  cop<- NULL
  ecop<-NULL
  qex<-c(quantile(EEx^2,qnt)-quantile(Ex^2,qnt),
         quantile(Ex,qnt),quantile(Ex^2,qnt))
  qey<-c(quantile(EEy^2,qnt)-quantile(Ey^2,qnt),
         quantile(Ey,qnt),quantile(Ey^2,qnt))
  
  
  
  if (copula){
    
    for (ix in quantile(x,qnt)){
      cx<-NULL
      for (iy in quantile(y,qnt)){
        cx<-c(cx,length(which(x<=ix & y <= iy))/length(x))
        
      }
      cop<-c(cop,(cx))
    }
    
    for (ix in quantile(x,qnt)){
      cx<-NULL
      for (iy in quantile(Ey,qnt)){
        cx<-c(cx,length(which(x<=ix & Ey <= iy))/length(x))
        
      }
      ecop<-c(ecop,(cx))
    }
    
  }
  
  lux<-length(unique(x))
  luy<-length(unique(y))
  
  N<-length(x)
  
  Ecdf.x=ecdf(x)(x) ## empirical cdf of D[,ef]
  Ecdf.y=ecdf(y)(y)
  
  asxy<-c(assoc(x,y),assoc(x,Ecdf.y),assoc(Ecdf.x,y),
          abs(cor(x2,y2))-abs(pcor1(x2,y2,Ex)),
          abs(cor(x1,y1))-abs(pcor1(x1,y1,Ey)))
  asey<-assoc(Ex,y2)
  asex<-assoc(Ex,x2)-asey
  aseyy<-assoc(Ey,y1)
  aseyx<-assoc(Ey,x1)-aseyy
  
  if (sd(Ex)>0.01)
    autx<-c(acf(Ex,plot=FALSE)$acf[2:3],pacf(Ex,plot=FALSE)$acf[1:2])
  else
    autx<-rep(1,4)
  
  if (sd(Ey)>0.01)
    auty<-c(acf(Ey,plot=FALSE)$acf[2:3],pacf(Ey,plot=FALSE)$acf[1:2])
  else
    auty<-rep(1,4)
  
  vxy<-varpred(x,y)
  vyx<-varpred(y,x)
  dv<- c(vxy,vxy-vyx)
  
  nX<-c(paste0("qx",1:length(qx)),paste0("qy",1:length(qy)),
        paste0("dv",1:length(dv)),paste0("cop",1:length(cop)),paste0("N",1:length(N)),
        paste0("mx",1:length(mx)),paste0("dx",1:length(dx)),
        paste0("my",1:length(my)),paste0("dy",1:length(dy)),
        paste0("qex",1:length(qex)),paste0("qey",1:length(qey)),
        paste0("lux",1:length(lux)),paste0("luy",1:length(luy)),
        paste0("asxy",1:length(asxy)),paste0("asex",1:length(asex)),
        paste0("asey",1:length(asey)),paste0("aseyx",1:length(aseyx)),
        paste0("aseyy",1:length(aseyy)),
        paste0("autx",1:length(autx)),paste0("auty",1:length(auty)))
  
  rX<-c(qx,qy,dv,cop,N, mx-my,dx-dy,my,dy,qex-qey,qey,lux,
        luy,asxy,asex,asey,aseyx, aseyy,autx,auty)
  
  names(rX)<-nX
  rX
  
}



varpred<-function(x,y,R=50){
  x<-scale(x)
  y<-scale(y)
  xh<-seq(-1.25,1.25,by=0.1)
  N<-length(x)
  P<-NULL
  beta<-NULL
  for (r in 1:R){
    #set.seed(r)
    ##    wr<-runif(N,0,5)
    ##    wr<-wr/sum(wr)
    Ir<-sample(1:N,round(4*N/5)) ##,replace=TRUE,prob=wr)
    xr<-x[Ir]
    yr<-y[Ir]
    
    px<-NULL
    for ( h in 1:length(xh)){
      sx<-sort(abs(xr-xh[h]),decreasing=FALSE,index=TRUE)$ix[1:min(10,length(xr))]
      px<-c(px,mean(yr[sx]))
      
      
    }
    
    P<-cbind(P,px)
    
  }
  
  
  mean(apply(P,1,sd))
  
  
}
