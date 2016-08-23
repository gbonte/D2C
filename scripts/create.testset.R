rm(list=ls())


noNodes<-c(20,100)
## range of number of nodes

N<-c(20,100)
## range of number of samples


NDAG.test=200

sdev<-c(0.2,1)

goParallel=FALSE
savefile<-TRUE
namefile<-paste("./data/testDAG",NDAG,type,"RData",sep=".")


testDAG<-new("simulatedDAG",NDAG=NDAG.test, N=N, noNodes=noNodes,
             functionType = c("linear","quadratic","sigmoid","kernel"),
             seed=101,sdn=sdev,quantize=c(FALSE),
             additive=c(TRUE,FALSE),goParallel=goParallel)


if (savefile)
  save(file=namefile,list=c("trainD2C","testDAG"))


