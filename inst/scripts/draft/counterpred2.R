rm(list=ls())
library(randomForest)
library(graph)
library(devtools)
#install_github("gbonte/D2C")
#install_github("gbonte/gbcode")
library(D2C)
library(gbcode)
load("../D2Cdata/trainD2Cp.1000.50.is.parent.RData")
trainD2CP<-trainD2C
load("../D2Cdata/trainD2Cp.1000.50.is.descendant.RData")
trainD2CD<-trainD2C
load("../D2Cdata/trainD2Cp.1000.50.is.ancestor.RData")

BACC<-NULL
WACC<-NULL
wACC<-NULL
ACC<-NULL
BNMSE<-NULL
WNMSE<-NULL
NMSE<-NULL
IACC<-NULL
ACC2<-NULL
NMSE2<-NULL
NMSE0<-NULL
ACC0<-NULL
for (r in 5:100000){

  set.seed(r+11)
  maxN=500
  n=sample(15:50,1)
  wgt = runif(n = 1,min = 0.85,max = 1)
  
  rDAG<-gRbase::random_dag(1:n,maxpar = 3,wgt)
  print("made rDAG")
  Edges= gRbase::edgeList(rDAG)
  
  Nodes<-nodes(rDAG)
  irDAG=igraph::graph.adjacency(as(rDAG,"matrix"))
  
  bestSP=0
  NX=NULL
  listN=NULL
  for (nn in Nodes){
    for (jj in setdiff(Nodes,nn)){
      SP=as.numeric(igraph::shortest.paths(irDAG, nn, jj, "out"))
      if ((!is.infinite(SP)) & SP >0){
        
        NX=nn
        NY=jj
        nascX=0
        for (nnn in Nodes){
          if (is.descendant(irDAG,nnn,NX)  ){
            nascX<-nascX+1     
          }}
          
        if (nascX>=5 &  SP>2)
          listN=c(listN,list(c(NX,NY,SP,nascX)))
      
    
  }}}
    
  
    if (length(listN)>=1 ){
      for (rr in 1:min(2,length(listN))){
        NX=listN[[rr]][1]
        NY=listN[[rr]][2]
        PaY<-NULL
        DescX<-NULL
        for (nn in Nodes){
          if (is.parent(irDAG,nn,NY) || is.parent(irDAG,NY,nn) ){
            PaY<-c(PaY,as.numeric(nn))      
          }
          if (is.descendant(irDAG,nn,NX)){
            DescX<-c(DescX,as.numeric(nn))      
          }
        }
        
        
        cat("r=", r, "NX=",NX,"NY=", NY,"LP=",listN[[rr]][3],
            "NascX=",listN[[rr]][4], "PaY=",PaY,"DescX=",DescX,"\n")
        NY<-as.numeric(NY) 
        NX=as.numeric(NX)
        H1 = function() return(H_Rn(1))
        H2 = function() return(H_Rn(2))
        Hs = function() return(H_sigmoid(1))
        Hk = function() return(H_kernel())
        N=0
        cnt=r
        while (N<50){
          set.seed(cnt)
          cnt=cnt+1
          HH=sample(c(Hk,Hs,H1),1)[[1]]
          DAG = new("DAG.network",network=rDAG,
        additive=FALSE,
          sdn=c(0.2,0.5),exosdn=runif(1,2,3),
          weights=c(1,2) ,H = HH)
        
        N=sample((maxN/2):maxN,1)
        
        X=compute(DAG,N)                              
        N=NROW(X)
        }
        delta=sample(X[,NX],N,rep=TRUE)
        #delta=sample(c(-1,1),rep=TRUE)+rnorm(2*N,sd=0.01)
        Xc0=counterfact(DAG,X,delta=delta,knocked=NX)
        print("out")
        
        Xc1=counterfact(DAG,X,delta=delta+sample(c(2,-2),N,rep=TRUE),knocked=NX)
      
        if (FALSE){
          graphics.off()
          colnames(X)<-1:n
          colnames(X)[NY]="Y"
          colnames(X)[NX]="X"
          dev.new()
          plot(irDAG)
          dev.new()
          par(mfrow=c(2,round(n/2)))
          for (i in setdiff(1:n,NY)){
            plot(X[,i],X[,NY],xlab=colnames(X)[i])
            points(Xc0[,i],Xc0[,NY],col="red")
            points(Xc1[,i],Xc1[,NY],col="green")
            
          }
          
          colnames(Xc0)<-1:n
          colnames(Xc0)[NY]="Y"
          colnames(Xc0)[NX]="X"
          #dev.new()
          #plot(data.frame(X))
          #dev.new()
          #plot(data.frame(Xc0),col="red")
         
          
        }
        if (FALSE){
        colnames(X)<-1:n
        colnames(X)[NY]="Y"
        colnames(X)[NX]="X"
        colnames(Xc0)<-1:n
        colnames(Xc0)[NY]="Y"
        colnames(Xc0)[NX]="X"
        XX<-rbind(X,Xc0)
        XX<-cbind(1:NROW(XX),XX)
        colnames(XX)[1]="ind"
        
        D=data.frame(XX)
        ggpairs(D,aes(color=ifelse(ind<N ,"green","black")))
        }
        
        N=NROW(X)
        print(N)
       
############### Compute the training observational error   
         

      Y0=Xc0[,NY]
      Y1=Xc1[,NY]
      
      cat(length(which(abs(Y1-Y0)<1e-3)),"sd=", sd(Y1-Y0),"\n")
      
      if (length(which(abs(Y1-Y0)<1e-3))>round(N/2))
        browser()
      

        
      } # for rr
    
  } # if ! is.null
}


