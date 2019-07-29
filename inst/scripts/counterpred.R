rm(list=ls())
library(graph)
library(D2C)
library(gbcode)

load("./data/trainD2C.1000.is.descendant.RData")
trainD2CD<-trainD2C
load("./data/trainD2C.500.is.ancestor.RData")


BACC<-NULL
WACC<-NULL
wACC<-NULL
ACC<-NULL
IACC<-NULL
ACC2<-NULL
for (r in 1:1000){
  set.seed(r+11)
  N=100
  n=sample(10:40,1)
  wgt = runif(n = 1,min = 0.95,max = 1)
  
  rDAG<-gRbase::random_dag(1:n,maxpar = 5,wgt) 
  Edges= gRbase::edgeList(rDAG)
  
  Nodes<-nodes(rDAG)
  topologicalOrder<-RBGL::tsort(rDAG)
  
  NY<-topologicalOrder[sample(1:n,1)] 
  
  irDAG=igraph::graph.adjacency(as(rDAG,"matrix"))
  plot(irDAG)
  
  
  NX<-NULL
  for (nn in topologicalOrder[round(n/2):(match(NY,topologicalOrder)-1)]){
    if (is.ancestor(irDAG,nn,NY)){
      NX<-as.character(nn)
      break
    }
  }
  if (! is.null(NX)){
    
    
    path<-topologicalOrder[match(NX,topologicalOrder):match(NY,topologicalOrder)]
    
    if (FALSE)
      for (nn in path){
        if (! is.ancestor(irDAG,nn,NY)){
          path=setdiff(path,nn)
        }
      }
    
    print(path)
    if (length(path)>1){
      NY<-as.numeric(NY) 
      NX=as.numeric(NX)
      
      
      
      
      DAG = new("DAG.network",network=rDAG,additive=TRUE,
                sdn=0.5,H = function() return(H_sigmoid(1)))
      
      X=compute(DAG,N)
      ancX=numeric(n)
      ancY=numeric(n)
      descX=numeric(n)
      for (i in 1:n){
        ancX[i] = predict(trainD2C,i,NX, X)$prob
        ancY[i] = predict(trainD2C,i,NY, X)$prob
        descX[i] = predict(trainD2CD,i,NX, X)$prob
        
       
      }
      sX=setdiff(sort(ancY-ancX+descX,decr=TRUE,index=TRUE)$ix[1:round(n/3)],c(which(ancY<0.5),NX,NY))
      
      
      Ipath=as.character(c(NX,sX,NY))
      cat("Ipath=",Ipath,"% inters=",length(intersect(sX,path))/length(path),"\n")
      delta=sample(c(-1,1),N,rep=TRUE)
      Xc=counterfact(DAG,X,delta=delta,knocked=NX)
      
      Xct=counterfact(DAG,X,delta=delta,knocked=NX)
      
      delta2=sample(c(-0.5,0.5),N,rep=TRUE)
      delta2=delta
      Xc2=counterfact(DAG,X,delta=delta2,knocked=NX)
      
      
      XX=X[,path]
      XXhat=XX*0
      XXhat[,1]=XX[,1]
      
      
      for (p in 2:NCOL(XX)){
        fs<-1:(p-1)
        if (length(fs)>5)
          fs<-union(1,fs[mrmr(XX[,fs],XX[,p])])
        ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],conPar=c(3,10),class=FALSE)
        XXhat[,p]=pred("rf",XX[,fs],XX[,p],XXhat[,fs],class=FALSE)
        e=XX[,p]-XXhat[,p]
        #cat(mean(e^2)/var(XX[,p]), " ")
      }
      e=XX[,NCOL(XX)]-XXhat[,NCOL(XX)]
      
      print(mean(e^2)/var(XX[,NCOL(XX)]))
      
      ###############@
      ###############
      ############@##### optimal CF
      XXct=Xct[,path]
      XXc=Xc[,path]
      XXhat=XX*0
      XXhat[,1]=XXc[,1]
      
      for (p in 2:(NCOL(XX))){
        fs<-1:(p-1)
        if (length(fs)>5)
          fs<-fs[mrmr(XX[,fs],XX[,p])]
        #fs<-1
        ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],linPar=c(5,10),cmbPar=10,class=FALSE)
        XXhat[,p]=pred("rf",XXct[,fs],XXct[,p],XXhat[,fs],class=FALSE)
        
      }
      
      acc=length(which((XXc[,NCOL(XXc)]-XX[,NCOL(XX)])*(XXhat[,NCOL(XXc)]-XX[,NCOL(XX)])>0))/N
      e=XXc[,NCOL(XXc)]-XXhat[,NCOL(XX)]
      BACC<-c(BACC,acc)
      cat("Best Acc=",mean(BACC), "NMSE=", mean(e^2)/var(XXc[,NCOL(XX)]), "\n")
      
      #########################################
      ############# ALLVAR
      Yc=Xc[,path[length(path)]]
      YY=XX[,path[length(path)]]
      XXhat=XX
      XXhat[,path[1]]=XXc[,path[1]]  ### change the column of the counterfactual
      
      fs<-1:NCOL(XX)
      if (length(fs)>5)
        fs<-fs[mrmr(XX[,fs],YY)]
      
      Yhat=pred("rf",XX[,fs],YY,XXhat[,fs],class=FALSE)
      
      acc=length(which((Yc-YY)*(Yhat-YY)>0))/N
      ## ITE= Yc-YY is the individualized treatment effect
      e=Yc-Yhat
      WACC<-c(WACC,acc)
      cat("Worst Acc=",mean(WACC), "NMSE=", mean(e^2)/var(Yc), "\n")
      
      #########################################
      ############# SELECTED CAUSES
      XXc=Xc[,path]
      XXhat=XX*0
      XXhat[,1]=XXc[,1]
      
      for (p in 2:(NCOL(XX))){
        fs<-1:(p-1)
        if (length(fs)>5)
          fs<-fs[mrmr(XX[,fs],XX[,p])]
        #fs<-1
        ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],linPar=c(5,10),cmbPar=10,class=FALSE)
        XXhat[,p]=pred("rf",XX[,fs],XX[,p],XXhat[,fs],class=FALSE)
        
      }
      
      acc=length(which((XXc[,NCOL(XXc)]-XX[,NCOL(XX)])*(XXhat[,NCOL(XXc)]-XX[,NCOL(XX)])>0))/N
      e=XXc[,NCOL(XXc)]-XXhat[,NCOL(XX)]
      ACC<-c(ACC,acc)
      cat("Acc=",mean(ACC), "NMSE=", mean(e^2)/var(XXc[,NCOL(XX)]), "\n")
      #########################################
      ############# INFERRED CAUSES
      XX=X[,Ipath]
      XXc=Xc[,Ipath]
      XXhat=XX*0
      XXhat[,1]=XXc[,1]
      
      for (p in 2:(NCOL(XX))){
        fs<-1:(p-1)
        if (length(fs)>5)
          fs<-fs[mrmr(XX[,fs],XX[,p])]
        #fs<-1
        ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],linPar=c(5,10),cmbPar=10,class=FALSE)
        XXhat[,p]=pred("rf",XX[,fs],XX[,p],XXhat[,fs],class=FALSE)
        
      }
      
      acc=length(which((XXc[,NCOL(XXc)]-XX[,NCOL(XX)])*(XXhat[,NCOL(XXc)]-XX[,NCOL(XX)])>0))/N
      e=XXc[,NCOL(XXc)]-XXhat[,NCOL(XX)]
      IACC<-c(IACC,acc)
      cat("Inferred Acc=",mean(IACC), "NMSE=", mean(e^2)/var(XXc[,NCOL(XX)]), "\n")
      #########################################
      if (FALSE){
      #########################################
      #########################################
      ############# SELECTED CAUSES + REWEIGHTED
      XX=X[,path]
      XXc=Xc[,path]
      XXhat=XX*0
      XXhat[,1]=XXc[,1]
      #D<-as.matrix(dist(cbind(XXc[,1],XX[,1])))
      
      for (i in 1:N){
        I<-sort(abs(XXc[i,1]-XX[,1]),decr=FALSE,index=TRUE)$ix[1:round(N/3)]
        for (p in 2:(NCOL(XX))){
          fs<-1:(p-1)
          if (length(fs)>5)
            fs<-fs[mrmr(XX[I,fs],XX[I,p])]
          #fs<-1
          ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],linPar=c(5,10),cmbPar=10,class=FALSE)
          XXhat[i,p]=pred("rf",XX[I,fs],XX[I,p],XXhat[i,fs],class=FALSE)
          
        }
      }
      acc=length(which((XXc[,NCOL(XXc)]-XX[,NCOL(XX)])*(XXhat[,NCOL(XXc)]-XX[,NCOL(XX)])>0))/N
      e=XXc[,NCOL(XXc)]-XXhat[,NCOL(XX)]
      wACC<-c(wACC,acc)
      cat("Weight Acc=",mean(wACC), "NMSE=", mean(e^2)/var(XXc[,NCOL(XX)]), "\n")
      }
      #########################################
      XX=X[,path]
      XXc=Xc2[,path]
      XXhat=XX*0
      XXhat[,1]=XXc[,1]
      
      for (p in 2:NCOL(XX)){
        fs<-1:(p-1)
        if (length(fs)>5)
          fs<-union(1,fs[mrmr(XX[,fs],XX[,p])])
        ##XXhat[,p]=pred("lazy",XX[,fs],XX[,p],XXhat[,fs],linPar=c(5,10),cmbPar=10,class=FALSE)
        XXhat[,p]=pred("rf",XX[,fs],XX[,p],XXhat[,fs],class=FALSE)
      }
      e=XXc[,NCOL(XXc)]-XXhat[,NCOL(XX)]
      acc=length(which((XXc[,NCOL(XXc)]-XX[,NCOL(XX)])*(XXhat[,NCOL(XXc)]-XX[,NCOL(XX)])>0))/N
      ACC2<-c(ACC2,acc)
      cat("Acc2=",mean(ACC2), "NMSE=", mean(e^2)/var(XXc[,NCOL(XX)]), "\n")
    }
  }
}


