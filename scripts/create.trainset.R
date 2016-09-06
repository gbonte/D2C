## Script that generates a training set for the D2C algorithm

rm(list=ls())


noNodes<-c(10,30)
## range of number of nodes

N<-c(50,200)
## range of number of samples

NDAG=200
## number of DAGs to be created and simulated

type="is.parent"

sdev<-c(0.2,0.5)
savefile<-TRUE
seed<-10
namefile<-paste("./data/trainD2C",NDAG,type,"RData",sep=".")

trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid","kernel"),
              seed=seed,sdn=sdev,quantize=c(TRUE,FALSE),maxpar.pc=c(0.05,0.3),
              additive=c(TRUE,FALSE),goParallel=FALSE)


descr.example<-new("D2C.descriptor",bivariate=TRUE,ns=3,acc=TRUE,lin=FALSE)

trainD2C<-new("D2C",sDAG=trainDAG,
              descr=descr.example,ratioEdges=0.5,
              max.features=10, type=type,goParallel=FALSE)



print(dim(trainD2C@X))
print(table(trainD2C@Y))
print(trainD2C@mod)


if (savefile)
  save(file=namefile,list=c("trainD2C","trainDAG","N","noNodes","seed"))

