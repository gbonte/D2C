## D2C: createTS.R
## D2C for time series
## The script
## 1) generates a training and a test set with different multivariate time series and
## 2) it learns a D2C classifier

rm(list=ls())
library(D2C)
library(doParallel)
type="is.parent"


set.seed(0)
noNodes<-6
## range of number of lags

N<-c(150,300)
## range of number of samples

NDAG=1000
## number of DAGs to be created and simulated
NDAG.test=500
nseries=c(5,100)
sdev<-c(0.1,0.3)

goParallel=TRUE
savefile<-TRUE
namefile<-"./data/traintestTSERIES.RData"
nfeat=30
maxs=50
ncores=1
if (goParallel){
  ncores=5
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
}

Sts=sample(1:25,10) 
Str=setdiff(1:25,Sts)


trainDAG<-new("simulatedTS",NDAG=NDAG, N=N, noNodes=noNodes,
              seed=10,sdn=sdev,goParallel=goParallel,nseries=nseries,
              typeser=Str)
cat("Computed trainDAG \n")
testDAG<-new("simulatedTS",NDAG=NDAG.test, N=N, noNodes=noNodes,
             seed=3101,sdn=sdev,
             goParallel=goParallel,typeser=Sts,
             nseries=nseries)

cat("Computed testDAG \n")


descr<-new("D2C.descriptor",bivariate=FALSE,ns=8,maxs=maxs,acc=TRUE,
           lin=TRUE,struct=TRUE, boot="rank",residual=FALSE,diff=FALSE,stabD=FALSE)

D2C<-new("D2C",sDAG=trainDAG,
         descr=descr,ratioEdges=0.15,
         max.features=nfeat, type=type,goParallel=goParallel,
         verbose=TRUE,npar=min(NDAG,ncores),rev=FALSE)


trainD2C<-makeModel(D2C,classifier="RF")


if (savefile)
    save(file=namefile,list=c("trainD2C","testDAG","noNodes"))
}




