rm(list=ls())
library(devtools)
install_github("gbonte/D2C")
require(gbcode)
library(D2C)
library(doParallel)
type="is.parent"

is.what<-function(iDAG,i,j){
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

set.seed(0)
noNodes<-c(15,100)

N<-c(50,100)
## range of number of samples

sdev<-c(0.1,0.5)
goParallel=TRUE

nfeat=20
maxs=10
ncores=1


NDAG=50000000 ## total number of DAGS
dN=5   ## number of DAGs per iteration
if (goParallel){
  ncores=10
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
  dN=3*ncores
}

  

iter=1
rep=1
while ( iter <=NDAG){
  
  set.seed(iter)

  if (rep%%2==0)
  trainDAG<-new("simulatedDAG",NDAG=dN, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"),
              seed=iter,sdn=sdev,verbose=TRUE,
              additive=c(TRUE,FALSE),goParallel=goParallel)
  else
  trainDAG<-new("simulatedTS",NDAG=dN, N=N+100, noNodes=c(5,10),
              seed=iter,sdn=sdev,verbose=TRUE,typeser=1:23,
              goParallel=goParallel,nseries=10)
	      
cat("Computed trainDAG \n")
  
  
  descr<-new("D2C.descriptor",bivariate=FALSE,ns=8,maxs=maxs,resid=TRUE,
             acc=TRUE,diff=FALSE,boot="rank",
             lin=TRUE,struct=FALSE,pq=c(0.1,0.5,0.9))

  
  D2C<-new("D2C",sDAG=trainDAG,
           descr=descr,ratioEdges=0.003,
           max.features=nfeat, type=type,goParallel=goParallel,
           verbose=TRUE,npar=min(NDAG,ncores))

   

  if (iter==1){
    allD2C<-D2C
  } else{
    allD2C@origX<-rbind(allD2C@origX,D2C@origX)
    allD2C@Y<-c(allD2C@Y,D2C@Y)
    cat("dim(allD2C@origX)=",dim(allD2C@origX),"NDAG computed=",iter+dN,"\n")
  }

  iter=iter+dN

    namefile<-paste("./data/trainD2Cp",length(allD2C@Y),max(noNodes),type,"RData",sep=".")
    trainD2C<-makeModel(allD2C,classifier="RF",EErep=5,max.features=50)
    save(file=namefile,list=c("trainD2C","allD2C","descr"))
    cat("SAVED ", namefile, "\n")
  

  rep=rep+1  
}


namefile<-paste("./data/trainD2Cp",NDAG,max(noNodes),type,"RData",sep=".")
trainD2C<-makeModel(allD2C,classifier="RF",EErep=10)
save(file=namefile,list=c("trainD2C","allD2C","descr"))
cat("SAVED ", namefile, "\n")




