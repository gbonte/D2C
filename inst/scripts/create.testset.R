rm(list=ls())
library(doParallel)
library(devtools)
install_github("gbonte/D2C")
library(D2C)

DAGSnames<-NULL

goParallel=TRUE
savefile<-TRUE
nDAGs<-100
namefile<-paste("./data/testDAGs",nDAGs,"RData",sep=".")
if (goParallel){
  ncores=20
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
  
}

for (ntests in 1:nDAGs){
  noNodes<-c(10,200)
  ## range of number of nodes
  
  N<-500
  ## range of number of samples
  
  NDAG.test=1000
  sdev<-c(0.2,1)
  
  if (ntests%%2==0)
    testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
                 functionType = c("linear","quadratic","sigmoid"),
                 seed=ntests,sdn=sdev,verbose=TRUE,
                 additive=c(TRUE,FALSE),goParallel=goParallel)
  else
    testDAG<-new("simulatedTS",NDAG=NDAG.test, N=N, noNodes=c(5,10),
                 seed=ntests,sdn=sdev,verbose=TRUE,typeser=1:23,
                 goParallel=goParallel,nseries=10)
  
  
  
  assign(paste("testDAG", ntests, sep = ""), testDAG)
  
  DAGSnames<-c(DAGSnames,paste("testDAG", ntests, sep = ""))
   
  if (savefile){
    save(file=namefile,list=c("nDAGs",DAGSnames))
  
    }
  
  
}