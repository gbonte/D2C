## Script that generates a training set for the D2C algorithm

rm(list=ls())


noNodes<-c(20,100)
## range of number of nodes

N<-c(20,100)
## range of number of samples

NDAG=50
## number of DAGs to be created and simulated

sdev<-c(0.2,1)
savefile<-TRUE
seed<-0
namefile<-paste("./data/trainD2C",NDAG,"RData",sep=".")

trainDAG<-new("simulatedDAG",NDAG=NDAG, N=N, noNodes=noNodes,
              functionType = c("linear","quadratic","sigmoid"),
              seed=seed,sdn=sdev,quantize=c(TRUE,FALSE),
              additive=c(FALSE),goParallel=FALSE)


descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=5,acc=TRUE,lin=FALSE)

trainD2C<-new("D2C",sDAG=trainDAG,
              descr=descr.example,ratioEdges=0.5,
              max.features=30, type=type,goParallel=FALSE)



print(dim(trainD2C@X))
print(table(trainD2C@Y))
print(trainD2C@mod)


if (savefile)
  save(file=namefile,list=c("trainD2C","trainDAG","N","noNodes","seed"))

