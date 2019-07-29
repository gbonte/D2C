rm(list=ls())
library(graph)
library(gbcode)
for (n in seq (15, 100, by=1)){
  N=1000
  wgt = runif(n = 1,min = 0.95,max = 1)
  
  rDAG<-gRbase::random_dag(1:n,maxpar = 5,wgt) 
  Edges= gRbase::edgeList(rDAG)
  
  Nodes<-nodes(rDAG)
  topologicalOrder<-RBGL::tsort(rDAG)
  ss<-sort(sample(1:n,2),decr=FALSE)
  NX<-topologicalOrder[ss[1]] 
  NY<-topologicalOrder[ss[2]] 
  
 
  NY<-as.numeric(NY) 
  NX=as.numeric(NX)
  print(NX)
  print(NY)
  irDAG=igraph::graph.adjacency(as(rDAG,"matrix"))
  plot(irDAG)
  
  
  
  inX<-NULL
  for (nn in setdiff(1:n,c(NX,NY))){
    if (is.parent(irDAG,as.character(nn),as.character(NX))
        && !is.parent(irDAG,as.character(nn),as.character(NY)) )
      inX<-c(inX,nn)
  }
  
  
  
  DAG = new("DAG.network",network=rDAG,additive=TRUE,sdn=0.5)
  
  torem=which(as.character(unlist(lapply(Edges,"[",2)))==NX)
  doDAG=rDAG
  if (length(torem)>0)
    for (r in torem)
      doDAG=graph::removeEdge(Edges[r][[1]][1],Edges[r][[1]][2],doDAG)
  
  idoDAG=igraph::graph.adjacency(as(doDAG,"matrix"))
  #plot(idoDAG)
  X=compute(DAG,N)
  X2=compute(DAG,N)
  
  Xc=counterfact(DAG,X,knocked=NX)
  browser()
  Yhat=pred("rf",X[,-NY],X[,NY],X2[,-NY],class=FALSE)
  Y=X2[,NY]
  
  cat("all tr:",mean((Yhat-Y)^2)/var(Y),"\n")
  
  Xk=compute(DAG,10*N,knocked=NX)
  #sX<-scale(rbind(X,Xk))
  
  #X<-sX[1:N,]
  #Xk<-sX[(N+1):NROW(sX),]
  Yhat=pred("rf",X[,-c(NY)],X[,NY],Xk[,-c(NY)],class=FALSE)
  Y=Xk[,NY]
  cat("all:",mean((Yhat-Y)^2)/var(Y),"\n")
  
  
  Yhat=pred("rf",X[,inX],X[,NY],Xk[,inX],class=FALSE)
  Y=Xk[,NY]
  cat("all+inX:",mean((Yhat-Y)^2)/var(Y),"\n")
  
  Yhat=pred("rf",X[,-c(inX,NY)],X[,NY],Xk[,-c(inX,NY)],class=FALSE)
  Y=Xk[,NY]
  cat("all-inX:",mean((Yhat-Y)^2)/var(Y),"\n")
  
  if (FALSE){
    fset=NX
    
    for (sz in 1:10){
      err=numeric(n-1)+Inf
      for (i in setdiff(1:n,c(fset,NY))){
        
        Yhat=pred("rf",X[,c(i,fset)],X[,NY],Xk[,c(i,fset)],class=FALSE)
        Y=Xk[,NY]
        err[i]=mean((Yhat-Y)^2)/var(Y)
        #cat("i=",i,igraph::shortest.paths(irDAG,as.character(NY),as.character(i),"in"),
        #   igraph::shortest.paths(idoDAG,as.character(NY),as.character(i),"in"),
        #  mean((Yhat-Y)^2)/var(Y),"\n")
      }
      fset<-c(fset,which.min(err))
      cat("fset=",fset, ":",min(err),"\n")
    }
    
    for (f in fset)
      cat("i=",f,igraph::shortest.paths(irDAG,as.character(NY),as.character(f),"in"),
          igraph::shortest.paths(idoDAG,as.character(NY),as.character(f),"in"),
          igraph::shortest.paths(irDAG,as.character(NX),as.character(f),"out"),
          igraph::shortest.paths(idoDAG,as.character(NX),as.character(f),"out"),
          "\n")
  }
}