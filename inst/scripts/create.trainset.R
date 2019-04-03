## Script that generates a training set for the D2C algorithm

rm(list=ls())
library(doParallel)
library(foreach)
library(D2C)
ncores=10
cl <- makeForkCluster(ncores, outfile='LOG.TXT')
registerDoParallel(cl)

noNodes<-c(10,100)
## range of number of nodes

N<-c(50,200)
## range of number of samples

NDAG=50
## number of DAGs to be created and simulated

type="is.descendant"

goParallel=TRUE
sdev<-c(0.2,0.5)
savefile<-TRUE
seed<-101
namefile<-paste("./data/trainD2C",NDAG,type,"RData",sep=".")

trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid","kernel"),
              seed=seed,sdn=sdev,quantize=c(TRUE,FALSE),maxpar.pc=c(0.05,0.3),
              additive=c(TRUE,FALSE),goParallel=goParallel)


descr.example<-new("D2C.descriptor",bivariate=TRUE,ns=3,acc=TRUE,lin=FALSE)

trainD2C<-new("D2C",sDAG=trainDAG,
              descr=descr.example,ratioEdges=0.6,
              max.features=15, type=type,goParallel=goParallel)



print(dim(trainD2C@X))
print(table(trainD2C@Y))
print(trainD2C@mod)


if (savefile)
  save(file=namefile,list=c("trainD2C","trainDAG","N","noNodes","seed"))

