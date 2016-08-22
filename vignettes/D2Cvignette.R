## ----echo=TRUE,results=TRUE----------------------------------------------
rm(list=ls())
library(D2C)
require(RBGL)
require(gRbase)



noNodes<-c(10,20)
## range of number of nodes

N<-c(50,200)
## range of number of samples

sd.noise<-c(0.2,1)
## range of values for standard deviation of additive noise 

NDAG=100
## number of DAGs to be created and simulated


trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"), 
              seed=0,sdn=sd.noise,quantize=c(TRUE,FALSE),verbose=FALSE)



## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@NDAG)

## ---- echo=TRUE----------------------------------------------------------
print(trainDAG@list.DAGs[[1]])
print(dim(trainDAG@list.observationsDAGs[[1]]))

