rm(list=ls())


noNodes<-c(40,120)
## range of number of nodes

N<-c(50,100)
## range of number of samples


NDAG.test=500

sdev<-c(0.2,1)

goParallel=FALSE
savefile<-TRUE
namefile<-paste("./data/testDAG",NDAG.test,"RData",sep=".")


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid","kernel"),
             seed=101,sdn=sdev,quantize=c(TRUE,FALSE),
             additive=c(TRUE,FALSE),goParallel=goParallel)


if (savefile)
  save(file=namefile,list=c("testDAG"))


