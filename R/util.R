## igraph
## g <- make_empty_graph(n = 5) %>%
##     add_edges(c(1,2, 3,2)) 


## This function allows the use of existing packages (e.g. pcalg) for the generation of
## random DAGs and related samples.
gendataDAG<-function(N,n,sdw=0.1){
  if (FALSE){
    D=data.frame(array(rnorm(N*n),c(N,n)))
  bn = gs(D)
  observedData<-rbn(bn,N,D)
  trueDAG<-as.graphNEL(bn)
  } else {
  trueDAG <- randomDAG(n, prob = runif(1,0.2,0.3), lB=0.1, uB=1)
  observedData <- scale(rmvDAG(N, trueDAG, 
                               errDist=sample(c("mix","cauchy","t4","mixt3"),1),mix=0.3))+array(rnorm(N*n,sd=sdw),c(N,n))
  }
  if (is.na(sum(observedData)))
    stop("error in gendataDAG")
  if (NCOL(observedData)!=length(nodes(trueDAG)))
    stop("error 2 in gendataDAG")
  list(DAG=trueDAG,data=observedData)
  
}

## high order correlation
HOC<-function(x,y,i,j){
  return(mean((x-mean(x))^i*(y-mean(y))^j)/((sd(x)^i)*(sd(y)^j)))
  
}

is.parent<-function(DAG,n1,n2){
  ## n1 : character name of node 1
  ## n2 : character name of node 2
  shortest.paths(DAG,n1,n2,"out")==1
  
}

is.child<-function(DAG,n1,n2){
  
  s=shortest.paths(DAG,n1,n2,"in")
  if (is.infinite(s))
    return(FALSE)
  return(s==1)
}

dagdistance<-function(DAG,n1,n2){
  
  dout=distances(DAG,n1,n2,mode="out")
  din=distances(DAG,n1,n2,mode="in")
  if (is.infinite(dout) & is.infinite(din) )
    return(sample(c(-5,5),1))
  if (is.infinite(dout) & is.finite(din) )
    return(-din)
  if (is.finite(dout) & is.infinite(din) )
    return(dout)
}


is.ancestor<-function(DAG,n1,n2){
  
  s=shortest.paths(DAG,n1,n2,"out")
  if (is.infinite(s))
    return(FALSE)
  return(s>=1)
}

is.descendant<-function(DAG,n1,n2){
  
  s=shortest.paths(DAG,n1,n2,"in")
  if (is.infinite(s))
    return(FALSE)
  return(s>=1)
  
  
}

is.mb<-function(DAG,n1,n2){
  
  is.child(DAG,n1,n2)||is.parent(DAG,n1,n2)
  
}


is.what<-function(iDAG,i,j,type){
  if (type=="is.mb")
    return(as.numeric(is.mb(iDAG,i,j)))
  if (type=="is.parent")
    return(as.numeric(is.parent(iDAG,i,j)))
  
  if (type=="is.child")
    return(as.numeric(is.child(iDAG,i,j)))
  if (type=="is.descendant")
    return(as.numeric(is.descendant(iDAG,i,j)))
  
  if (type=="is.ancestor")
    return(as.numeric(is.ancestor(iDAG,i,j)))
  
}

Sigm<-function(x,W=1){
  E=exp(x)
  return(2*W*E/(1+E)-W)
}

rankrho<-function(X,Y,nmax=5,regr=FALSE,first=NULL){
  ## mutual information ranking
  ## 17/10/11
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  if (var(Y)<0.01)
    return(1:nmax)
  X<-scale(X)
  
  Iy<-numeric(n)
  if (!regr){
    Iy<-cor2I2(corXY(X,Y))
  } else {
    for (i in 1:n)
      Iy[i]<-abs(regrlin(X[,i],Y)$beta.hat[2])
  }
  
  if (m>1)
    Iy<-apply(Iy,1,mean)
  
  
  return(sort(c(Iy), decreasing=TRUE, index.return=TRUE)$ix[1:nmax])
  
  
}




quantization<-function(x,nbin=1){
  if (nbin==1)
    return(as.numeric(cut(x, breaks = c(min(x)-1,median(x),max(x)+1))))
  else
    return(as.numeric(cut(x, breaks = c(min(x)-1,quantile(x,c(0.25,0.5,0.75)),max(x)+1))))
  
}

H_sigmoid <- function(n=2)
{
  a = runif(n+1,min = -1,max = 1)
  f <- function(x)
  {
    X  = x^(0:n)
    return ( -1+2/(1+exp(mean(X * a))))
  }
  return(Vectorize(f))
}



H_Rn <- function(n){
  a = runif(n+1,min = -1,max = 1)
  f <- function(x)
  {
    X  = x^(0:n)
    return ( sum(X * a))
  }
  return(Vectorize(f))
}


kernel.fct<- function(X,knl=anovadot(sigma=runif(1,0.5,2),
                                     degree=sample(1:2,1)),lambda=0.1){
  save.seed <- get(".Random.seed", .GlobalEnv)
  N<-NROW(X)
  set.seed(as.numeric(Sys.time()))
  Y<-rnorm(N,sd=1)
  K<-kernelMatrix(knl,X)
  Yhat<-K%*%ginv(K+lambda*N*diag(N))%*%Y
  assign(".Random.seed", save.seed, .GlobalEnv)
  return(Yhat)
}

H_kernel <- function()
{
  
  f <- function(x)
  {
    
    return (kernel.fct(x))
  }
  return(Vectorize(f))
}


pcor1<-function(x,y,z){
  ## partial correlation cor(x,y|z)
  
  
  if (is.numeric(z)){
    rho.xy<-cor(x,y,"pairwise.complete.obs")
    rho.xz<-cor(x,z,"pairwise.complete.obs")
    rho.yz<-cor(z,y,"pairwise.complete.obs")
    if (is.na(rho.xz+rho.yz+rho.xy))
      return(0)
    if (rho.xz==1 | rho.yz==1)
      return(0)
    rho<-(rho.xy-rho.xz*rho.yz)/(sqrt(1-min(rho.xz^2,0.99))*sqrt(1-min(rho.yz^2,0.99)))
    return(rho)
  } else {
    stop("z should be numeric")
    
  }
  
  
  
}


corDC<-function(X,Y){
  ## correlation continuous matrix and discrete vector
  ## NB: the notion of sign has no meaning in this case. Mean of absolute values is taken
  ## 14/11/2011
  
  if (!is.factor(Y))
    stop("This is not the right function. Y is not a factor !!")
  
  N<-NROW(X)
  L<-levels(Y)
  
  if( length(L)==2)
    lL<-1
  else
    lL<-length(L)
  
  cxy<-NULL
  for (i in 1:lL){
    yy<-numeric(N)
    ind1<-which(Y==L[i])
    ind2<-setdiff(1:N,ind1)
    yy[ind1]<-1
    cxy<-cbind(cxy,abs(cor(X,yy)))
  }
  
  apply(cxy,1,mean)
}


Icond<-function(x,y=NULL,z,lambda=0){
  ## conditional  information cor(x,y|z)
  
  
  
  ## numeric z
  if (is.numeric(z)){
    if (is.vector(x))
      return(cor2I2(pcor1(x,y,z)))
    X<-x
    n<-NCOL(X)
    Ic<-array(0,c(n,n))
    for (i in 1:(n-1))
      for (j in (i+1):n){
        Ic[i,j]<-Icond(X[,i],X[,j],z)
        Ic[j,i]<-Ic[i,j]
      }
    return(Ic)
    
  }
  ## factor z and vectors x and y
  if (! is.null(y)){
    L<-levels(z)
    lL<-length(L)
    w<-numeric(lL)
    for (i in 1:lL)
      w[i]<-length(which(z==L[i]))
    w<-w/sum(w)
    
    Ic<-NULL
    for (i in 1:lL){
      ind1<-which(z==L[i])
      Ic<-c(Ic,cor2I2(cor(x[ind1],y[ind1])))
    }
    
    return(as.numeric(w*Ic))
  }
  
  ## factor z and matrix x
  X<-x
  n<-NCOL(X)
  L<-levels(z)
  lL<-length(L)
  w<-numeric(lL)
  for (i in 1:lL)
    w[i]<-length(which(z==L[i]))
  w<-w/sum(w)
  
  Ic<-array(0,c(n,n))
  W<-0
  for (i in 1:lL){
    ind1<-which(z==L[i])
    
    
    if (length(ind1)>8){
      
      
      Ic<-Ic+w[i]*cor2I2(cor.shrink(X[ind1,],lambda=lambda,verbose=F))
      W<-W+w[i]
    }
  }
  
  
  return(Ic/W)
  
}









ppears<-function(r.hat,N,S=0){
  n<-length(r.hat)
  p<-numeric(n)
  
  for (i in 1:n){
    z<-abs(0.5*(log(1+r.hat[i])-log(1-r.hat[i])))*sqrt(N[i]-S-3)
    
    p[i]<-pnorm(z,lower.tail=F)
    
    
  }
  p
}

corXY<-function(X,Y){
  ## correlation continuous matrix and continuous/discrete vectormatrix
  
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  cXY<-array(NA,c(n,m))
  
  for (i in 1:m){
    if (m==1)
      YY<-Y
    else
      YY<-Y[,i]
    if (is.numeric(YY)){
      cXY[,i]<-cor(X,YY,use="pairwise.complete.obs")
    } else {
      cXY[,i]<-corDC(X,YY)
    }
  }
  cXY
}


cor2I2<-function(rho){
  rho<-pmin(rho,1-1e-5)
  rho<-pmax(rho,-1+1e-5)
  -1/2*log(1-rho^2)
  
  
}


lazy.pred<- function(X,Y,X.ts,class=FALSE,return.more=FALSE,
                     conPar=3,linPar=5,cmbPar=10){
  
  n<-NCOL(X)
  N<-NROW(X)
  
  if (class){ ## classification
    l.Y<-levels(Y)
    L<-length(l.Y)
    u<-unique(Y)
    
    if (length(u)==1){
      P<-array(0,c(NROW(X.ts),L))
      colnames(P)<-l.Y
      P[,u]<-1
      out.hat<-factor(rep(as.character(u),length(X.ts)),levels=l.Y)
      return(list(pred=out.hat,prob=P))
    }
    
    if (L==2) {
      
      
      stop("not supported")
    } else {
      algo="lazy"
      
      stop("not supported")
      
    }
  } else { ## regression
    d<-data.frame(cbind(Y,X))
    names(d)[1]<-"Y"
    
    names(d)[2:(n+1)]<-paste("x",1:n,sep="")
    
    
    
    mod<-lazy(Y~.,d,control=lazy.control(distance="euclidean",
                                         conIdPar=conPar,
                                         linIdPar=linPar,
                                         cmbPar=cmbPar))
    if (is.vector(X.ts) & n>1)
      X.ts<-array(X.ts,c(1,n))
    d.ts<-data.frame(X.ts)
    
    names(d.ts)<-names(d)[2:(n+1)]
    
    if (!return.more){
      ll<- predict(mod,d.ts)
      return(ll$h)
      
    } else {
      ll<- predict(mod,d.ts,S.out=TRUE,k.out=FALSE)
      return(ll)
    }
  }
  
}



rf.pred<- function(X,Y,X.ts,class=FALSE,...){
  
  n<-NCOL(X)
  N<-NROW(X)
  if (is.vector(X.ts) & n>1){
    N.ts<-1
    X.ts<-array(X.ts,c(1,n))
  }  else {
    if (n==1)
      N.ts<-length(X.ts)
    else
      N.ts<-nrow(X.ts)
  }
  
  if (n==1){
    X<-array(X,c(N,1))
    X.ts<-array(X.ts,c(N.ts,1))
  }
  #set.seed(N*n)
  d<-data.frame(Y,X)
  names(d)[1]<-"Y"
  
  mod.rf<-randomForest(Y~.,data=d,...)
  d.ts<-data.frame(X.ts)
  names(d.ts)[1:(n)]<-names(d)[2:(n+1)]
  p<-predict(mod.rf,d.ts,type="response")
  
  if (class){
    P<-predict(mod.rf,d.ts,type="prob")
    return(list(pred=p,prob=P))
  } else {
    
    
    return(p)
  }
  
}


#' mIMR (minimum Interaction max Relevance) filter
#' @param X :  input matrix
#' @param Y : output vector
#' @param nmax : number of returned features
#' @param init : if TRUE it makes a search in the space of pairs of features to initialize the ranking, otherwise the first ranked feature is the one with the highest mutual information with the output
#' @param lambda : weight \eqn{0 \le \lambda \le 1} of the interaction term
#' @param spouse.removal : TRUE OR FALSE. if TRUE it removes the spouses before ranking
#' @param caus :   if \code{caus =1} it prioritizes causes otherwise (\code{caus=-1}) it prioritizes effects
#' @return ranked vector of \code{nmax} indices of features
#' @description Filter  based on information theory which aims to prioritise direct causal relationships in feature selection problems where the ratio between the number of features and the number of samples is high. The approach is based on the notion of interaction which is informative about the relevance of an input subset as well as its causal relationship with the target.
#' @examples
#'  set.seed(0)
#' N<-500
#' n<-5
#' X<-array(rnorm(N*n),c(N,n))
#' Y<-X[,1]-3*X[,3]+4*X[,2]+rnorm(N,sd=0.5)
#' Z1<-Y+rnorm(N,sd=0.5)
#' ## effect 1
#' Z2<-2*Y+rnorm(N,sd=0.5)
#' ## effect 2
#' most.probable.causes<-mimr(cbind(X,Z1,Z2),Y,nmax=3,init=TRUE,spouse=FALSE,lambda=1)
#' ## causes are in the first three columns of the feature dataset
#' most.probable.effects<-mimr(cbind(X,Z1,Z2),Y,nmax=3,init=TRUE,spouse=FALSE,lambda=1,caus=-1)
#' ## effects are in the last two columns of the feature dataset
#' @references Bontempi G., Meyer P.E. (2010) Causal filter selection in microarray data. ICML10
#' @export
mimr<-function(X,Y,nmax=5,
               init=FALSE,lambda=0.5,
               spouse.removal=TRUE,
               caus=1){
  if (var(Y)<0.01)
    return(1:nmax)
  NMAX<-nmax
  m<-NCOL(Y) # number of outputs
  n<-NCOL(X)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY<-var(Y)
  CY<-corXY(X,Y)
  Iy<-cor2I2(CY)
  subset<-1:n
  pv.rho<-ppears(c(CY),N+numeric(n))
  if (spouse.removal){
    pv<-ppears(c(CY),N+numeric(n))
    s<-sort(pv,decreasing=TRUE,index.return=TRUE)
    hw<-max(1,min(n-nmax,length(which(s$x>0.05))))
    spouse<-s$ix[1:hw]
    subset<-setdiff(1:n,s$ix[1:hw])
    X<-X[,subset]
    Iy<-Iy[subset]
    n<-NCOL(X)
  }
  
  
  CCx<-cor(X)
  Ix<-cor2I2(CCx)
  ## mutual information
  Ixx<-Icond(X,z=Y,lambda=0.02)
  ## conditional information
  Inter<-array(NA,c(n,n))
  
  if (init){
    max.kj<--Inf
    for (kk in 1:(n-1)){
      for (jj in (kk+1):n){
        Inter[kk,jj]<- (1-lambda)*(Iy[kk]+Iy[jj])+caus*lambda*(Ixx[kk,jj]-Ix[kk,jj])
        Inter[jj,kk]<-Inter[kk,jj]
        if (Inter[kk,jj]>max.kj){
          max.kj<-Inter[kk,jj]
          subs<-c(kk,jj)
        }
      }
    }
  } else {
    subs<-which.max(Iy)
  }
  
  if (nmax>length(subs)){
    last.subs<-0
    for (j in length(subs):min(n-1,NMAX-1)){
      mrmr<-numeric(n)-Inf
      if (length(subs)<(n-1)){
        if (length(subs)>1){
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*apply(-Ix[subs,-subs]+Ixx[subs,-subs],2,mean)
        } else {
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*(-Ix[subs,-subs]+Ixx[subs,-subs])
        }
      } else {
        mrmr[-subs]<-Inf
      }
      s<-which.max(mrmr)
      subs<-c(subs,s)
    }
  }
  
  ra<-subset[subs]
  
  if (nmax>length(ra))
    ra<-c(ra,setdiff(1:orign,ra))
  
  ra
  
}


mimr2<-function(X,Y1,Y2,nmax=5,
                init=FALSE,lambda=0.5,
                spouse.removal=TRUE,
                caus=1){
  if (sd(Y1)<1e-5 || sd(Y2)<1e-5)
    return(list(fs1=1:nmax, fs2=1:nmax))
  NMAX<-nmax
  m<-NCOL(Y1) # number of outputs
  n<-NCOL(X)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY1<-var(Y1)
  CY1<-corXY(X,Y1)
  Iy1<-cor2I2(CY1)
  HY2<-var(Y2)
  CY2<-corXY(X,Y2)
  Iy2<-cor2I2(CY2)
  subset1<-1:n
  pv.rho1<-ppears(c(CY1),N+numeric(n))
  if (spouse.removal){
    pv1<-ppears(c(CY1),N+numeric(n))
    s1<-sort(pv1,decreasing=TRUE,index.return=TRUE)
    hw<-max(1,min(n-nmax,length(which(s1$x>0.05))))
    spouse<-s1$ix[1:hw]
    subset1<-setdiff(1:n,s1$ix[1:hw])
    
    pv2<-ppears(c(CY2),N+numeric(n))
    s2<-sort(pv2,decreasing=TRUE,index.return=TRUE)
    hw<-max(1,min(n-nmax,length(which(s2$x>0.05))))
    spouse<-s2$ix[1:hw]
    subset2<-setdiff(1:n,s2$ix[1:hw])
    
    subset<-union(subset1,subset2)
    
    X<-X[,subset]
    Iy1<-Iy1[subset]
    Iy2<-Iy2[subset]
    n<-NCOL(X)
  }
  
  
  CCx<-cor(X)
  Ix<-cor2I2(CCx)
  ## mutual information
  Ixx1<-Icond(X,z=Y1,lambda=0.02)
  Ixx2<-Icond(X,z=Y2,lambda=0.02)
  ## conditional information
  Inter1<-array(NA,c(n,n))
  Inter2<-array(NA,c(n,n))
  if (init){
    max.kj<--Inf
    for (kk in 1:(n-1)){
      for (jj in (kk+1):n){
        Inter1[kk,jj]<- (1-lambda)*(Iy1[kk]+Iy1[jj])+caus*lambda*(Ixx1[kk,jj]-Ix[kk,jj])
        Inter1[jj,kk]<-Inter1[kk,jj]
        if (Inter1[kk,jj]>max.kj){
          max.kj<-Inter1[kk,jj]
          subs1<-c(kk,jj)
        }
      }
    }
    max.kj<--Inf
    for (kk in 1:(n-1)){
      for (jj in (kk+1):n){
        Inter2[kk,jj]<- (1-lambda)*(Iy2[kk]+Iy2[jj])+caus*lambda*(Ixx2[kk,jj]-Ix[kk,jj])
        Inter2[jj,kk]<-Inter2[kk,jj]
        if (Inter2[kk,jj]>max.kj){
          max.kj<-Inter2[kk,jj]
          subs2<-c(kk,jj)
        }
      }
    }
    
  } else {
    subs1<-which.max(Iy1)
    subs2<-which.max(Iy2)
  }
  
  if (nmax>length(subs1)){
    
    for (j in length(subs1):min(n-1,NMAX-1)){
      mrmr1<-numeric(n)-Inf
      mrmr2<-numeric(n)-Inf
      if (length(subs1)<(n-1)){
        if (length(subs1)>1){
          mrmr1[-subs1]<- (1-lambda)*Iy1[-subs1]+caus*lambda*apply(-Ix[subs1,-subs1]-Ix[subs2,-subs1]+Ixx1[subs1,-subs1],2,mean)
          mrmr2[-subs2]<- (1-lambda)*Iy2[-subs2]+caus*lambda*apply(-Ix[subs2,-subs2]-Ix[subs1,-subs2]+Ixx2[subs2,-subs2],2,mean)
        } else {
          mrmr1[-subs1]<- (1-lambda)*Iy1[-subs1]+caus*lambda*(-Ix[subs1,-subs1]-Ix[subs2,-subs1]+Ixx1[subs1,-subs1])
          mrmr2[-subs2]<- (1-lambda)*Iy2[-subs2]+caus*lambda*(-Ix[subs2,-subs2]-Ix[subs1,-subs2]+Ixx2[subs2,-subs2])
        }
      } else {
        mrmr1[-subs1]<-Inf
        mrmr2[-subs2]<-Inf
      }
      s1<-which.max(mrmr1)
      subs1<-c(subs1,s1)
      s2<-which.max(mrmr2)
      subs2<-c(subs2,s2)
    }
  }
  
  ra1<-subset[subs1]
  ra2<-subset[subs2]
  
  if (nmax>length(ra1)){
    ra1<-c(ra1,setdiff(1:orign,ra1))
    ra2<-c(ra2,setdiff(1:orign,ra2))
  }
  if (any(is.na(ra1))|| any(is.na(ra2))){
    stop("error in mimr2")
  }
  list(fs1=ra1,fs2=ra2)
  
}


mimrmat<-function(X,Y,nmax=5,
                  init=TRUE,lambda=0.5,
                  spouse.removal=TRUE,R=100){
  if (var(Y)<0.01)
    return(1:nmax)
  
  m<-NCOL(Y) # number of outputs
  n<-NCOL(X)
  nmax<-min(n-2,nmax)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY<-var(Y)
  CY<-corXY(X,Y)
  Iy<-cor2I2(CY)
  subset<-1:n
  pv.rho<-ppears(c(CY),N+numeric(n))
  if (spouse.removal){
    pv<-ppears(c(CY),N+numeric(n))
    s<-sort(pv,decreasing=TRUE,index.return=TRUE)
    hw<-max(1,min(n-nmax,length(which(s$x>0.05))))
    spouse<-s$ix[1:hw]
    subset<-setdiff(1:n,s$ix[1:hw])
    X<-X[,subset]
    Iy<-Iy[subset]
    n<-NCOL(X)
  }
  CCx<-cor(X)
  Ix<-cor2I2(CCx)
  ## mutual information
  Ixx<-Icond(X,z=Y,lambda=0.02)
  f<-function(x,y){
    if (x==1 && y==1)
      return(1)
    return(-1)
  }
  best=-Inf
  for (r in 1:R){
    if (runif(1)<0.5 || r <10)
      C<-sample(c(-1,1),n,replace=TRUE)
    else{
      C=bestC
      i=sample(n,1)
      C[i]=-C[i]
    }
    wC=which(C>0)
    Sy=sum(Iy[wC])
    MC=outer(C,C,Vectorize(f))
    S=MC*Ixx-MC*Ix
    diag(S)=0
    
    if ((sum(S)+Sy)>best){
      best=sum(S)+Sy
      bestC=C
      print(best)
    }
    
  }
  subs=sort(bestC,decr=TRUE,index=TRUE)$ix[1:nmax]
  
  ra<-subset[subs]
  
  if (nmax>length(ra))
    ra<-c(ra,setdiff(1:orign,ra))
  
  ra
}

mimreff<-function(X,Y,nmax=5,
                  init=TRUE,lambda=0.5,
                  spouse.removal=TRUE){
  if (var(Y)<0.01)
    return(1:nmax)
  
  m<-NCOL(Y) # number of outputs
  n<-NCOL(X)
  nmax<-min(n-2,nmax)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY<-var(Y)
  CY<-corXY(X,Y)
  Iy<-cor2I2(CY)
  subset<-1:n
  pv.rho<-ppears(c(CY),N+numeric(n))
  if (spouse.removal){
    pv<-ppears(c(CY),N+numeric(n))
    s<-sort(pv,decreasing=TRUE,index.return=TRUE)
    hw<-max(1,min(n-nmax,length(which(s$x>0.05))))
    spouse<-s$ix[1:hw]
    subset<-setdiff(1:n,s$ix[1:hw])
    X<-X[,subset]
    Iy<-Iy[subset]
    n<-NCOL(X)
  }
  
  
  CCx<-cor(X)
  Ix<-cor2I2(CCx)
  ## mutual information
  Ixx<-Icond(X,z=Y,lambda=0.02)
  ## conditional information
  Inter<-array(NA,c(n,n))
  
  max.kj<--Inf
  for (kk in 1:(n-1)){
    for (jj in (kk+1):n){
      Inter[kk,jj]<- (1-lambda)*(Iy[kk]+Iy[jj])+lambda*(Ixx[kk,jj]-Ix[kk,jj])
      Inter[jj,kk]<-Inter[kk,jj]
      if (Inter[kk,jj]>max.kj){
        max.kj<-Inter[kk,jj]
        causes<-c(kk,jj)
      }
    }
  }
  
  
  subs<-NULL
  
  while (length(subs)< nmax){
    score<-numeric(n)-Inf
    subs2<-c(subs,causes)
    if (length(subs2)>=(n-1))
      score[-subs2]<- (1-lambda)*Iy[-subs2]-lambda*mean(-Ix[subs2,-subs2]+Ixx[subs2,-subs2])
    else
      score[-subs2]<- (1-lambda)*Iy[-subs2]-lambda*apply(-Ix[subs2,-subs2]+Ixx[subs2,-subs2],2,mean)
    
    s<-which.max(score)
    subs<-c(subs,s)
  }
  ra<-subset[subs]
  if (nmax>length(ra))
    ra<-c(ra,setdiff(1:orign,ra))
  
  ra
}




mrmr<-function(X,Y,nmax=5,first=NULL,all=FALSE,back=FALSE,lambda=1,categ=FALSE){
  ## mRMR filter
  # 17/10/11
  
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  if (categ && is.factor(Y)){
    Iy<-numeric(n)
    Ix<-array(0,c(n,n))
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        Ix[i,j]<-mean(c(mutinf(factor(X[,i]),factor(X[,j])),
                        mutinf(factor(X[,j]),factor(X[,i]))))
      }
      Iy[i]<-mutinf(factor(X[,i]),Y)
    }
  }else {
    
    X<-scale(X)
    Iy<-cor2I2(corXY(X,Y))
    
    CCx<-cor(X,use="pairwise.complete.obs")
    Ix<-cor2I2(CCx)
    
  }
  
  subs<-which.max(Iy)
  for (j in length(subs):min(n-1,nmax)){
    mrmr<-numeric(n)-Inf
    if (length(subs)<(n-1)){
      if (length(subs)>1){
        mrmr[-subs]<- Iy[-subs]+lambda*apply(-Ix[subs,-subs],2,mean)
      } else {
        
        mrmr[-subs]<- Iy[-subs]+lambda*(-Ix[subs,-subs])
        
      }
    } else {
      mrmr[-subs]<-Inf
    }
    
    s<-which.max(mrmr)
    sortmrmr<-sort(mrmr,decreas=TRUE,index=TRUE)$ix[1:(n-length(subs))]
    allfs<-c(subs,sortmrmr)
    subs<-c(subs,s)
    
  }
  
  
  if (back){  ## backward reordering based on linear regression
    nsubs<-NULL
    while (length(subs)>1){
      pd<-numeric(length(subs))
      
      for (ii in 1:length(subs))
        pd[ii]<-regrlin(X[,setdiff(subs,subs[ii])],Y)$MSE.emp
      
      nsubs<-c(subs[which.min(pd)],nsubs)
      subs<-setdiff(subs,subs[which.min(pd)])
      
      
    }
    subs<-c(subs,nsubs)
  }
  
  
  if (all){
    return(allfs)
  } else {
    return(subs[1:nmax])
  }
  
}



mrmr2<-function(X,Y1,Y2,nmax=5,first=NULL,all=FALSE,lambda=1,categ=FALSE){
  ## mRMR2 filter
  # 26/11/19
  
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y1)
  
  if (categ && is.factor(Y1)){
    Iy1<-numeric(n)
    Ix<-array(0,c(n,n))
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        Ix[i,j]<-mean(c(mutinf(factor(X[,i]),factor(X[,j])),
                        mutinf(factor(X[,j]),factor(X[,i]))))
      }
      Iy1[i]<-mutinf(factor(X[,i]),Y1)
      Iy2[i]<-mutinf(factor(X[,i]),Y2)
    }
    
    
  }else {
    
    X<-scale(X)
    Iy1<-cor2I2(corXY(X,Y1))
    Iy2<-cor2I2(corXY(X,Y2))
    CCx<-cor(X,use="pairwise.complete.obs")
    Ix<-cor2I2(CCx)
    
  }
  
  subs1<-which.max(Iy1)
  subs2<-which.max(Iy2)
  if (subs2==subs1)
    subs2=sort(Iy2,decr=TRUE,index=TRUE)$ix[2]
  for (j in length(subs1):min(n-1,nmax)){
    mrmr1<-numeric(n)-Inf
    mrmr2<-numeric(n)-Inf
    if (length(subs1)<(n-1)){
      if (length(subs1)>1){
        mrmr1[-subs1]<- Iy1[-subs1]+lambda*apply(-Ix[subs1,-subs1],2,mean)+lambda*apply(-Ix[subs2,-subs1],2,mean)
        mrmr2[-subs2]<- Iy2[-subs2]+lambda*apply(-Ix[subs2,-subs2],2,mean)+lambda*apply(-Ix[subs1,-subs2],2,mean)
      } else {
        
        mrmr1[-subs1]<- Iy1[-subs1]+lambda*(-Ix[subs1,-subs1])+lambda*(-Ix[subs2,-subs1])
        mrmr2[-subs2]<- Iy2[-subs2]+lambda*(-Ix[subs1,-subs2])+lambda*(-Ix[subs1,-subs2])
        
      }
    } else {
      mrmr1[-subs1]<-Inf
      mrmr2[-subs2]<-Inf
    }
    
    s1<-which.max(mrmr1)
    s2<-which.max(mrmr2)
    sortmrmr1<-sort(mrmr1,decreas=TRUE,index=TRUE)$ix[1:(n-length(subs1))]
    sortmrmr2<-sort(mrmr2,decreas=TRUE,index=TRUE)$ix[1:(n-length(subs2))]
    allfs1<-c(subs1,sortmrmr1)
    subs1<-c(subs1,s1)
    
    allfs2<-c(subs2,sortmrmr2)
    subs2<-c(subs2,s2)
    
  }
  
  
  if (all){
    return(list(allfs1,allfs2))
  } else {
    return(list(fs1=subs1[1:nmax],fs2=subs2[1:nmax]))
  }
  
}


assoc <-function(x,y){
  c(abs(cor(x,y)),cor.test(x,y)$p.value)
  
}


#' Balanced Error Rate
#' @param Ytrue :  binary numeric vector (made of 0 or 1) of real classes
#' @param Yhat : binary numeric vector (made of 0 or 1) of predicted classes
#' @description The balanced error rate is the average of the errors on each class: BER = 0.5*(FP/(TN+FP) + FN/(FN+TP)).
#' @return Balanced Error Rate \eqn{0 \le } BER \eqn{ \le 1}
#' @export
BER<-function(Ytrue,Yhat){
  
  if (!(is.numeric(Ytrue) & is.numeric(Yhat)))
    stop("BER accepts only numeric values")
  TN<-length(which(Yhat==0 & Ytrue==0))
  FN<-length(which(Yhat==0 & Ytrue==1))
  TP<-length(which(Yhat==1 & Ytrue==1))
  FP<-length(which(Yhat==1 & Ytrue==0))
  
  
  b1<-FP/(TN+FP)
  b2<-FN/(FN+TP)
  if (is.na(b1))
    b1<-0
  if (is.na(b2))
    b2<-0
  return(0.5*(b1+b2))
  
  
  
}


#' AUC
#' @author Gianluca Bontempi  \email{gbonte@@ulb.ac.be}
#' @references Handbook \emph{Statistical foundations of machine learning} available in \url{http://www.ulb.ac.be/di/map/gbonte/mod_stoch/syl.pdf}
#' @description AUC
#' @details AUC
#' @title AUC
#' @name AUC
#' @export
#'
#' @param  y: real value
#' @param  yhat: predicted probability
#' @return AUC
#' @export
#' @examples
#' ## random prediction
#' AUC(round(runif(100)),rnorm(100))
#'
AUC<-function(y,yhat){
  
  p<-prediction(yhat,y)
  p<-performance(p,"auc")
  
  mean(unlist(slot(p,"y.values")),na.rm=T)
}



#### MakeEmbedded ####
#' Embed a multivariate time series in input output form
#' @author Gianluca Bontempi  \email{gbonte@@ulb.ac.be}
#' @references \url{mlg.ulb.ac.be}
#' @title Embedding of multivariate time series
#' @param ts: multivariate time series [no. observations,no. variates]
#' @param n [no.var]: vector of embedding orders
#' @param delay [no.var]: vector of delays
#' @param hor [no.var]: vector of predicted horizons (hor=1 boils down to one-step-ahed prediction)
#' @param w: index of variables appearing in $out
#'
#' @return list with
#' \itemize{
#' \item{inp}: embedded inputs of all variables (n[i] columns for series ts[i])
#' \item{out}: outputs of variates whose index is in w}
#' @examples
#' TS<-array(1:25,c(5,5))
#' ## embeds a 5 variate time series with embedding orders equal to 2 for all variables for one-step ahead prediction
#' MakeEmbedded(TS,n=rep(2,5),delay=rep(0,5),hor=rep(1,5),w=1:5)
#'
#'
MakeEmbedded<-function(ts, n, delay,hor=1,w=1){
  
  no.data<-NROW(ts)
  no.var<-NCOL(ts)
  a<-NROW(n)
  b<-NCOL(n)
  if (a!=no.var)
    stop('Error in the size of embedding n')
  if (length(delay)!=no.var)
    stop('Error in the size of delay')
  if (length(hor)!=length(w))
    stop('Error in the size of horizon hor')
  N<-no.data-max(n)-max(delay)
  
  Input<-array(0,c(N,sum(n)))
  Output<-array(0,c(N,sum(hor)))
  
  for (i in 1:N) {
    for (j in 1:no.var) {
      k<-1:n[j]
      Input[i,sum(n[1:j-1])+k]<-ts[i+n[j]-k+max(n)-n[j]+max(delay)-delay[j],j]
      
      for (ww in 1:length(w)){
        if (ww==1)
          iw<-0
        else
          iw<-sum(hor[1:(ww-1)])
        
        Output[i,(iw+1):(sum(hor[1:ww]))]<-numeric(hor[ww])+NA
        M<-min(no.data,(i+max(n)+max(delay)+hor[ww]-1))
        
        
        Output[i,(iw+1):(iw+M-(i+max(n)+max(delay))+1)]<-ts[(i+max(n)+max(delay)):M,w[ww]]
        
      }
    }
    
  }
  
  list(inp=Input,out=Output)
}



genTS<-function(nn,NN,sd=0.5,num=1){
  
  
  n=4  ## max embedding order 
  Y=rnorm(nn,sd=0.1)
  ep=0
  th0=rnorm(1)
  fs<-sample(0:(n),4)
  state=0
  if (num>0){
    for (i in 1:NN){
      N=length(Y)
      
      if (num==1){
        nfs=2
        e=rnorm(1)
        Y=c(Y,-0.4*(3-Y[N-fs[1]]^2)/(1+Y[N-fs[1]]^2)+0.6*(3-(Y[N-fs[2]]-0.5)^3)/(1+(Y[N-fs[2]]-0.5)^4)+sd*(e+th0*ep))
        ep=e
      }
      if (num==2){
        nfs=2
        e=rnorm(1)
        Y=c(Y,(0.4-2*exp(-50*Y[N-fs[1]]^2))*Y[N-fs[1]]+(0.5-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sd*(e+th0*ep))
        ep=e
      }
      if (num==3){
        nfs=3
        e=rnorm(1)
        Y=c(Y,1.5 *sin(pi/2*Y[N-fs[1]])-sin(pi/2*Y[N-fs[2]])+sin(pi/2*Y[N-fs[3]])+sd*(e+th0*ep))
        ep=e
      }
      if (num==4){
        nfs=2
        e=rnorm(1)
        Y=c(Y,2*exp(-0.1*Y[N-fs[1]]^2)*Y[N-fs[1]]-exp(-0.1*Y[N-fs[2]]^2)*Y[N-fs[2]]+sd*(e+th0*ep))
        ep=e
      }
      if (num==5){
        nfs=1
        Y=c(Y,-2*Y[N-fs[1]]*max(0,sign(-Y[N-fs[1]]))+0.4*Y[N-fs[1]]*max(0,sign(Y[N-fs[1]]))+sd*rnorm(1))
      }
      if (num==6){
        nfs=2
        Y=c(Y,0.8*log(1+3*Y[N-fs[1]]^2)-0.6*log(1+3*Y[N-fs[2]]^2)+sd*rnorm(1))
      }
      if (num==7){
        nfs=2
        Y=c(Y,1.5 *sin(pi/2*Y[N-fs[1]])-sin(pi/2*Y[N-fs[2]])+sd*rnorm(1))
      }
      if (num==8){
        nfs=2
        Y=c(Y,(0.5-1.1*exp(-50*Y[N-fs[1]]^2))*Y[N-fs[1]]+(0.3-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sd*rnorm(1))
      }
      if (num==9){
        nfs=2
        Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1))
      }
      if (num==10){
        nfs=1
        Y=c(Y,sign(Y[N-fs[1]])+sd*rnorm(1))
      }
      if (num==11){
        nfs=1
        Y=c(Y,0.8*Y[N-fs[1]]-0.8*Y[N-fs[1]]/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1))
      }
      
      if (num==12){
        nfs=2
        Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sd*rnorm(1))
      }
      if (num==13){
        nfs=1
        fs=0
        Y=c(Y,min(1,max(0,3.8*Y[N-fs[1]]*(1-Y[N-fs[1]])+sd*rnorm(1))))
      }
      if (num==14){
        nfs=2
        fs=0:1
        ## Henon map
        Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sd*rnorm(1))
      }
      if (num==15){
        nfs=1
        ## TAR map
        if (Y[N-fs[1]]<1){
          Y=c(Y,-0.5*Y[N-fs[1]]+sd*rnorm(1))
        } else {
          Y=c(Y,0.4*Y[N-fs[1]]+sd*rnorm(1))
        }
      }
      if (num==16){
        nfs=1
        ## TAR2 map
        if (state==1){
          Y=c(Y,-0.5*Y[N-fs[1]]+sd*rnorm(1))
        } else {
          Y=c(Y,0.4*Y[N-fs[1]]+sd*rnorm(1))
        }
        if (runif(1)>0.9)
          state=1-state
      }
      if (num==17){
        nfs=4
        #GARCH 
        Y=c(Y,sqrt(0.000019+0.846*((Y[N-fs[1]])^2+0.3*(Y[N-fs[2]])^2+0.2*(Y[N-fs[3]])^2+0.1*(Y[N-fs[4]])^2) )*sd*rnorm(1))
        
      }
      
      if (any(is.nan(Y) | abs(Y)>1000))
        stop("error")
      
    } ## for i
    
  }
  
  if (num<0){
    fs<- 1:(-num)
    nfs=length(fs)
    ord=-num
    Cf=rnorm(ord)
    ma=rnorm(2)
    #min(Mod(polyroot(c(1, -model$ar))))
    while (any(Mod(polyroot(c(1,-Cf)))<=1))
      Cf=rnorm(ord)
    Y<-c(arima.sim(n = NN, list(ar = Cf, ma = ma,sd = sd)))  
  }
  
  fs=fs[1:nfs]
  
  Y=scale(Y[nn:length(Y)])
  if (any(is.nan(Y) ))
    stop("error")
  M=MakeEmbedded(array(Y,c(length(Y),1)),n=nn,delay=0,hor=rep(1,1),w=1:1)
  netwDAG<-new("graphNEL", nodes=as.character(1:nn), edgemode="directed")
  
  for (j in 1:(nn-max(fs)-1)){
    for (f in fs){
      
      netwDAG <- addEdge(as.character(j+f+1), as.character(j), netwDAG, 1)
    }
  }
  
  list(D=M$inp,DAG=netwDAG)
}

genSTAR<-function(n, nn,NN,sdev=0.5,num=1,loc=2,verbose=FALSE){
  # n: number of time series
  # nn: max lag
  ## NN: number of samples
  ## loc : size neghborhood
  ##
  ## Returns:
  ## D [NN,n*nn] observation dataset
  ## DAG associated DAG
  ## nodes 1:nn= series 1
  ## nodes nn+1:2*nn= series 2
  if (nn<5)
    stop("Too few lags")
  
  ep=0
  doNeigh=list()
  for (i in 1:n){
    if (runif(1)<2/n){
      doNeigh[[i]]=sample(-1:1,sample(1:2,1))
      if (i==1)
        doNeigh[[i]]=sample(0:1,sample(1:2,1))
      if (i==n)
        doNeigh[[i]]=sample(-1:0,sample(1:2,1))
    } else {
      doNeigh[[i]]=0
    }
  }
  eold=numeric(n)
  eold2=numeric(n)
  th0=rnorm(1)
  fs<-1:4 #sample(0:(nn-2),4)
  ## subset of lags 
  state=0
  print(num)
  Y=array(rnorm(n*nn,sd=0.1),c(nn,n))
  rep=0
  repeat{
    
    for (ii in 1:NN){
      N=NROW(Y)
      y=numeric(n)
      if (num==1){
        nfs=2
        for (i in 1:n){ ## number of time series
          neigh=i+doNeigh[[i]] ## to which series the antecedents belong 
          
          e=rnorm(1)
          y[i]=-0.4*(3-mean(Y[N-fs[1],neigh])^2)/(1+mean(Y[N-fs[1],neigh])^2)+
            0.6*(3-(mean(Y[N-fs[2],neigh])-0.5)^3)/(1+(mean(Y[N-fs[2],neigh])-0.5)^4)+sdev*e
          
        }
        
      }
      if (num==2){
        nfs=2
        
        #Y=c(Y,(0.4-2*exp(-50*Y[N-fs[1]]^2))*Y[N-fs[1]]+(0.5-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sdev*(e+th0*ep))
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          
          e=rnorm(1)
          y[i]=(0.4-2*exp(-50*mean(Y[N-fs[1],neigh])^2))*mean(Y[N-fs[1],neigh])+
            (0.5-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sdev*e
          
        }
      }
      if (num==3){
        nfs=3
        #Y=c(Y,1.5 *sin(pi/2*Y[N-fs[1]])-sin(pi/2*Y[N-fs[2]])+sin(pi/2*Y[N-fs[3]])+sdev*(e+th0*ep))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          
          e=rnorm(1)
          y[i]=1.5 *sin(pi/2*mean(Y[N-fs[1],neigh]))-
            sin(pi/2*mean(Y[N-fs[2],neigh]))+
            sin(pi/2*mean(Y[N-fs[3],neigh]))+sdev*e
          
        }
      }
      
      if (num==4){
        nfs=2
        #  Y=c(Y,2*exp(-0.1*Y[N-fs[1]]^2)*Y[N-fs[1]]-exp(-0.1*Y[N-fs[2]]^2)*Y[N-fs[2]]+sdev*(e+th0*ep))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=2*exp(-0.1*mean(Y[N-fs[1],neigh])^2)*mean(Y[N-fs[1],neigh])-
            exp(-0.1*mean(Y[N-fs[2],neigh])^2)*mean(Y[N-fs[2],neigh])+sdev*e
          
        }
      }
      
      if (num==5){
        nfs=1
        #Y=c(Y,-2*Y[N-fs[1]]*max(0,sign(-Y[N-fs[1]]))+0.4*Y[N-fs[1]]*max(0,sign(Y[N-fs[1]]))+sdev*rnorm(1))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=-2*mean(Y[N-fs[1],neigh])*max(0,sign(-mean(Y[N-fs[1],neigh])))+
            0.4*mean(Y[N-fs[1],neigh])*max(0,sign(mean(Y[N-fs[1],neigh])))+sdev*e
          
        }
      }
      if (num==6){
        nfs=2
        #Y=c(Y,0.8*log(1+3*Y[N-fs[1]]^2)-0.6*log(1+3*Y[N-fs[2]]^2)+sdev*rnorm(1))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=0.8*log(1+3*mean(Y[N-fs[1],neigh])^2)-0.6*log(1+3*mean(Y[N-fs[2],neigh])^2)+sdev*e
          
        }
      }
      if (num==7){
        nfs=2
        #y[i]=(0.4-2*cos(40*mean(Y[N-5,neigh]))*exp(-30*mean(Y[N-5,neigh])^2))*mean(Y[N-5,neigh])+(0.5-0.5*exp(-50*mean(Y[N-9,neigh])^2))*mean(Y[N-9,neigh])
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=(0.4-2*cos(40*mean(Y[N-fs[1],neigh]))*exp(-30*mean(Y[N-fs[1],neigh])^2))*mean(Y[N-fs[1],neigh])+
            (0.5-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sdev*e
          
          
        }
      }
      if (num==8){
        nfs=2
        #Y=c(Y,(0.5-1.1*exp(-50*Y[N-fs[1]]^2))*Y[N]+(0.3-0.5*exp(-50*Y[N-fs[2]]^2))*Y[N-fs[2]]+sdev*rnorm(1))
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=0.5-1.1*exp(-50*mean(Y[N-fs[1],neigh])^2)*mean(Y[N-fs[1],neigh])+
            (0.3-0.5*exp(-50*mean(Y[N-fs[2],neigh])^2))*mean(Y[N-fs[2],neigh])+sdev*e
          
        }
        
      }
      
      if (num==9){
        nfs=2
        ##Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sdev*rnorm(1))
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=0.3*mean(Y[N-fs[1],neigh])+0.6*mean(Y[N-fs[2],neigh])+
            (0.1-0.9*mean(Y[N-fs[1],neigh])+
               0.8*mean(Y[N-fs[2],neigh]))/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sdev*e
        }
        
      }
      if (num==10){
        nfs=1
        ##Y=c(Y,sign(Y[N-fs[1]])+sdev*rnorm(1))
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=sign(mean(Y[N-fs[1],neigh]))+sdev*e
          
        }
      }
      if (num==11){
        nfs=1
        ##Y=c(Y,0.8*Y[N-fs[1]]-0.8*Y[N-fs[1]]/(1+exp(-10*Y[N-fs[1]]))+sdev*rnorm(1)) 
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=0.8*mean(Y[N-fs[1],neigh])-0.8*mean(Y[N-fs[1],neigh])/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sdev*e
          
        }
      }
      if (num==12){
        nfs=2
        ##  Y=c(Y,0.3*Y[N-fs[1]]+0.6*Y[N-fs[2]]+(0.1-0.9*Y[N-fs[1]]+0.8*Y[N-fs[2]])/(1+exp(-10*Y[N-fs[1]]))+sdev*rnorm(1))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=0.3*mean(Y[N-fs[1],neigh])+0.6*mean(Y[N-fs[2],neigh])+
            (0.1-0.9*mean(Y[N-fs[1],neigh])+0.8*mean(Y[N-fs[2],neigh]))/(1+exp(-10*mean(Y[N-fs[1],neigh])))+sdev*e
        }
        
      }
      if (num==13){
        nfs=1
        ##  Y=c(Y,min(1,max(0,3.8*Y[N-fs[1]]*(1-Y[N-fs[1]])+sdev*rnorm(1))))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=min(1,max(0,3.8*mean(Y[N-fs[1],neigh])*(1-mean(Y[N-fs[1],neigh]))+sdev*e))
          
        }
      }
      if (num==14){
        nfs=2
        ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sdev*rnorm(1))
        
        ## AR
        if (NROW(Y)<=nn){
          W=numeric(max(fs[1:nfs])+1)
          
          W[fs[1:nfs]+1]=rnorm(nfs)
          while (any(Mod(polyroot(c(1,-W)))<=1))
            W[fs[1:nfs]+1]=rnorm(nfs)
        }
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]= W[fs[1]+1]*mean(Y[N-fs[1],neigh]) + W[fs[2]+1]*mean(Y[N-fs[2],neigh])+sdev*e
          
        }
        
      }
      if (num==15){
        nfs=1
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          if (mean(Y[N-fs[1],neigh])<1){
            y[i]=-0.5*mean(Y[N-fs[1],neigh])+sdev*e
          } else {
            y[i]=0.4*mean(Y[N-fs[1],neigh])+sdev*e
          }
        }
      }
      
      if (num==16){
        nfs=1
        ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sdev*rnorm(1))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          if (abs(mean(Y[N-fs[1],neigh]))<=1){
            y[i]=0.9*mean(Y[N-fs[1],neigh])+sdev*e
          } else {
            y[i]=-0.3*mean(Y[N-fs[1],neigh])+sdev*e
          }
        }
      } 
      if (num==17){
        nfs=1
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          if (state==1){
            y[i]=-0.5*mean(Y[N-fs[1],neigh])+sdev*rnorm(1)
          } else {
            y[i]=0.4*mean(Y[N-fs[1],neigh])+sdev*rnorm(1)
          }
          if (runif(1)>0.9)
            state=1-state
        }
      }
      
      if (num==18){
        
        nfs=4
        #GARCH 
        ##Y=c(Y,sqrt(0.000019+0.846*((Y[N-fs[1]])^2+0.3*(Y[N-fs[2]])^2+0.2*(Y[N-fs[3]])^2+0.1*(Y[N-fs[4]])^2) )*sdev*rnorm(1))
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]=sqrt(0.000019+0.846*((mean(Y[N-fs[1],neigh]))^2+
                                      0.3*(mean(Y[N-fs[2],neigh]))^2+
                                      0.2*(mean(Y[N-fs[3],neigh]))^2+
                                      0.1*(mean(Y[N-fs[4],neigh]))^2) )*sdev*e
          
          
        }
      }
      
      if (num==19){
        nfs=1
        e=numeric(n)
        # y=0.7 y(t-1)*e(t-2)+e(t)
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e[i]=rnorm(1)
          y[i]=0.7*(mean(Y[N-fs[1],neigh]))*sdev*eold2[i]+sdev*e[i]
        }
        eold2=eold
        eold=e
      }
      
      if (num==20){
        nfs=2
        e=numeric(n)
        # y=0.4 y(t-1)-0.3 y(t-2)+0.5*y(t-1)*e(t-1)+e(t)
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e[i]=rnorm(1)
          y[i]=0.4*(mean(Y[N-fs[1],neigh]))-0.3*(mean(Y[N-fs[2],neigh]))+
            0.5*(mean(Y[N-fs[1],neigh]))*sdev*eold[i]+sdev*e[i]
        }
        eold2=eold
        eold=e
      }
      
      if (num==21){
        nfs=1
        e=numeric(n)
        # y=0.7 abs(y(t-1))/(abs(y(t-1))+2)+e
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e[i]=rnorm(1)
          y[i]=0.7*abs(mean(Y[N-fs[1],neigh]))/ (abs(mean(Y[N-fs[1],neigh]))+2)+sdev*e[i]
        }
        eold2=eold
        eold=e
      }
      
      if (num==22){
        nfs=1
        e=numeric(n)
        # y=0.9 y(t-1)- 0.8 y(t-1)/(1+exp(-10 y(t-1))) +e
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e[i]=rnorm(1)
          y[i]=0.9*mean(Y[N-fs[1],neigh])-
            0.8*mean(Y[N-fs[1],neigh])/ (1+exp(-10*mean(Y[N-fs[1],neigh])))+sdev*e[i]
          
        }
        eold2=eold
        eold=e
      }
      if (num==23){
        nfs=2
        e=numeric(n)
        # y=0.3 y(t-1)+ 0.6 y(t-2)+ (0.1-0.9 y(t-1)+0.8 y(t-2))/(1+exp(-10 y(t-1))) +e
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e[i]=rnorm(1)
          y[i]=0.3*mean(Y[N-fs[1],neigh]) +0.6 *mean(Y[N-fs[2],neigh])+(0.1-0.9*mean(Y[N-fs[1],neigh])+
                                                                          0.8*mean(Y[N-fs[2],neigh]))/ (1+exp(-10*mean(Y[N-fs[1],neigh])))+sdev*e[i]
        }
        eold2=eold
        eold=e
      }
      
      if (num==24){
        nfs=3
        ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sdev*rnorm(1))
        
        ## AR
        if (NROW(Y)<=nn){
          W=numeric(max(fs[1:nfs])+1)
          
          W[fs[1:nfs]+1]=rnorm(nfs)
          while (any(Mod(polyroot(c(1,-W)))<=1))
            W[fs[1:nfs]+1]=rnorm(nfs,sd=0.5)
        }
        
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]= W[fs[1]+1]*mean(Y[N-fs[1],neigh]) + W[fs[2]+1]*mean(Y[N-fs[2],neigh])+ 
            W[fs[3]+1]*mean(Y[N-fs[3],neigh])+sdev*e
          
        }
        
      }
      
      if (num==25){
        nfs=1
        ##   Y=c(Y,1-1.4*Y[N-fs[1]]*Y[N-fs[1]] + 0.3*(Y[N-fs[2]])+0.001*sdev*rnorm(1))
        
        ## AR
        if (NROW(Y)<=nn){
          W=numeric(max(fs[1:nfs])+1)
          
          W[fs[1:nfs]+1]=rnorm(nfs)
          while (any(Mod(polyroot(c(1,-W)))<=1))
            W[fs[1:nfs]+1]=rnorm(nfs)
        }
        
        for (i in 1:n){
          neigh=i+doNeigh[[i]]
          e=rnorm(1)
          y[i]= W[fs[1]+1]*mean(Y[N-fs[1],neigh]) +sdev*e
          
        }
        
      }
      
      Y<-rbind(Y,y)
      
      
    } ## for i
    if (! (any(is.nan(Y) | abs(Y)>10000))){
      break
    } else  {
      rep=rep+1
      Y=array(rnorm(n*nn,sd=0.1/rep),c(nn,n))
      if (rep>5)
        num<-sample(1:5,1)
     
    }
   
    
  }## repeat
  
  
  
  
  fs=fs[1:nfs]
  if (any(apply(Y,2,sd)<1e-3*sdev))
    stop("genSTAR error: too small sd in Y ")
  Y=scale(Y[nn:NROW(Y),])
  if (any(is.nan(Y) ))
    stop("genSTAR error: nan Y")
  M=MakeEmbedded(Y,n=numeric(n)+nn,delay=numeric(n),hor=rep(1,n),w=1:n)
  netwDAG<-new("graphNEL", nodes=as.character(1:(n*nn)), edgemode="directed")
  if (NCOL(M$inp)!=(n*nn))
    stop("genstar NCOL(M$inp)!=(n*nn)")
  fs=sort(fs)
  for (i in 1:n){ ## for each series
    for (j in 1:nn){ ## for each node in the series
      for (f in fs){ ## for each parent
        if ((j+f)<nn){
          Neigh=i+doNeigh[[i]] ## it defines the series where the parent is 
          
          for (neigh in Neigh){ ## for each parent
            netwDAG <- addEdge(as.character((neigh-1)*nn+j+f+1), as.character((i-1)*nn+j), netwDAG, 1)
            if (verbose){
             
              cat((neigh-1)*nn+j+f+1,"->",(i-1)*nn+j,"\n")
            }
          }
        }
      }
    }
  }
  
  if (any(apply(M$inp,2,sd)<(sdev*1e-3)))
    stop("genSTAR error: too small sd in M$inp ")
  list(D=M$inp,DAG=netwDAG,fs=fs,Y=Y,doNeigh=doNeigh)
}

## return the set of temporal possible causes for the effect
timecauses<-function(nD,nn,ind.effect){
  #nn=6 
  n=nD/nn
  #ind.effn2mfrow()ect=10
  tt<-ind.effect%%nn
  nn
  ind.causes<-NULL
  
  for (i in 1:n){
    if (tt<nn)
      ind.causes<-c(ind.causes,(i-1)*nn+seq(tt+1,nn))
  }
  return(ind.causes)
}
