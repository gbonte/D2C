rm(list=ls())
library(doParallel)
library(devtools)
install_github("gbonte/D2C")
library(D2C)
noNodes<-c(10,100)
## range of number of nodes

N<-c(50,100)
## range of number of samples


NDAG.test=10

sdev<-c(0.2,1)

goParallel=FALSE
savefile<-TRUE
namefile<-paste("./data/testDAG",NDAG.test,max(noNodes),"RData",sep=".")

if (goParallel){
  ncores=20
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
  
}

testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid","kernel"),
             seed=101,sdn=sdev,quantize=c(TRUE,FALSE),
             additive=c(TRUE,FALSE),goParallel=goParallel)


if (savefile)
  save(file=namefile,list=c("testDAG"))

print("SAVED")

