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
BACC2<-NULL
WACC<-NULL
WACC2<-NULL
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
for (r in 1:100000){
  
  set.seed(as.numeric(Sys.time())+r)
  maxN=500
  n=sample(10:30,1)
  wgt = runif(n = 1,min = 0.85,max = 1)
  
  rDAG<-gRbase::random_dag(1:n,maxpar = 5,wgt)
  
  nodes(rDAG)<-as.character(1:n)
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
        
        if (  SP>2 & nascX>3)
          listN=c(listN,list(c(NX,NY,SP,nascX)))
        
        
      }}}
  
  
  if (length(listN)>=1 ){
    for (rr in 1:min(5,length(listN))){
      NX=listN[[rr]][1]
      NY=listN[[rr]][2]
      PaY<-NULL
      ChY<-NULL
      DescX<-NULL
      for (nn in Nodes){
        if (is.parent(irDAG,nn,NY)  ){
          PaY<-c(PaY,as.numeric(nn))      
        }
        if ( is.parent(irDAG,NY,nn) ){
          ChY<-c(ChY,as.numeric(nn))      
        }
        if (is.descendant(irDAG,nn,NX)){
          DescX<-c(DescX,as.numeric(nn))      
        }
      }
      
      
      cat("r=", r, "NX=",NX,"NY=", NY,"LP=",listN[[rr]][3],
          "NascX=",listN[[rr]][4], "PaY=",PaY, "ChY=",ChY,"DescX=",DescX,"\n")
      NY<-as.numeric(NY) 
      NX=as.numeric(NX)
      H1 = function() return(H_Rn(1))
      H2 = function() return(H_Rn(2))
      Hs = function() return(H_sigmoid(1))
      Hk = function() return(H_kernel())
      N=0
      cnt=0
      L=0
      N=maxN
      while (L<100){
        set.seed(as.numeric(Sys.time())+cnt)
        cnt=cnt+1
        HH=sample(c(Hs,H1,H2),1)[[1]]
        DAG = new("DAG.network",network=rDAG,
                  additive=sample(c(TRUE,FALSE),1),
                  sdn=c(0.1,0.2),exosdn=runif(1,0.5,2),
                  weights=c(0.5,2) ,H = HH)
        
        N=maxN
        
        X=compute(DAG,N)                              
        N=NROW(X)
        if (N>(maxN/2)){
         if ( all(apply(X,2,sd)>0.08)) {
          set.seed(as.numeric(Sys.time())+cnt)
          DX=runif(1,0.25,0.5)
          
          #delta=sample(c(-1,1),rep=TRUE)+rnorm(2*N,sd=0.01)
          Xc0=counterfact(DAG,X,delta=sample(X[,NX],N,rep=TRUE),knocked=NX)
          
          delta=runif(N,DX,2*DX)*sample(c(-1,1),N,rep=TRUE)
          Xc1=counterfact(DAG,X,delta=delta ,knocked=NX)
          Y0=Xc0[,NY]
          Y1=Xc1[,NY]
          DY=unique(round(Y1-Y0,2))
          if (length(DY)<10)
            L=0
          else
            L=length(which(abs(DY)>DX/3))*(max(abs(DY))<4)
          if (L>100)
            cat("L=",L,"DX=",DX,"sd=",sd(Y1-Y0),summary(DY),"\n")
         
        } else {
          N=maxN
        }}
        else {
          N=maxN
        }
        if (cnt>10){
          N=0
          print("BREAK")
          break
        }
      }
      
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
      
      
      if (N>50){
        X=scale(X)
        Xc0=scale(Xc0,center=attr(X,'scaled:center'),scale=attr(X,'scaled:scale'))
        Xc1=scale(Xc1,center=attr(X,'scaled:center'),scale=attr(X,'scaled:scale'))
        Y0=Xc0[,NY]
        Y1=Xc1[,NY]
        
        
        ancX=numeric(n)
        parY=numeric(n)
        descX=numeric(n)
        for (i in 1:n){
          if (i!=NX)
            ancX[i] = predict(trainD2CP,i,NX, X)$prob
          if (i !=NY)
            parY[i] = predict(trainD2CP,i,NY, X)$prob
          if (i !=NX)
            descX[i] = predict(trainD2CD,i,NX, X)$prob
          
          
        }
        sX=union(setdiff(sort(parY*(1-descX),decr=TRUE,index=TRUE)$ix[1],c(NX,NY)),
                 setdiff(which(parY>0.5),c(NX,NY)))
        
        ## set of variables parent of Y and no descendants of X
        px=which.max(ancX)
        if (max(ancX)>0.5)
          px=sort(ancX,decr=TRUE,index=TRUE)$ix[1:2]
        ## set of variables parent of X 
        cat("sX=",sX,"parY=",which(parY>0.5),"\n")
        
        
        ############### Compute the training observational error   
        
        #  if (length(sX)>0)
        #    Yhat=pred("rf",X[,c(sX,NX)],X[,NY],X[,c(sX,NX)],class=FALSE)
        #  else
        #    Yhat=pred("rf",X[,NX],X[,NY],X[,NX],class=FALSE)
        #  eps=X[,NY]-Yhat
        
        
        
       
        
        XX0=X
        XX0[,NX]=Xc0[,NX]
        
        XX1=X
        XX1[,NX]=Xc1[,NX]
        ############### counterfactual0
        
        Yhat0=pred("rf",X[,NX],X[,NY],XX0[,NX],class=FALSE)
        Yhat1=pred("rf",X[,NX],X[,NY],XX1[,NX],class=FALSE)
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))/abs(DD)
        
        ACC0<-c(ACC0,acc)
        
        NMSE0<-c(NMSE0, mean(abs(e)) )
        cat("Cnt Acc0=",mean(ACC0), "NMSE=",mean(NMSE0), "\n")
        ############### counterfactual
        
        if (length(sX)>0){
          Yhat0=pred("rf",X[,c(sX,NX)],X[,NY],XX0[,c(sX,NX)],class=FALSE)
          Yhat1=pred("rf",X[,c(sX,NX)],X[,NY],XX1[,c(sX,NX)],class=FALSE)
        }
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))/abs(DD)
        ACC<-c(ACC,acc)
        NMSE<-c(NMSE, mean(abs(e)))
        cat("Cnt Acc=",mean(ACC), "NMSE=",mean(NMSE), "\n")
        ############### 
        if (length(px)>0){
          
          Xbck=X[,px]
          if (length(sX)>0) {
            Yhat0=pred("rf",cbind(X[,c(sX,NX)],Xbck),X[,NY],cbind(XX0[,c(sX,NX)],Xbck),class=FALSE)
            Yhat1=pred("rf",cbind(X[,c(sX,NX)],Xbck),X[,NY],cbind(XX1[,c(sX,NX)],Xbck),class=FALSE)
          }else{
            Yhat0=pred("rf",cbind(X[,c(NX)],Xbck),X[,NY],cbind(XX0[,c(NX)],Xbck),class=FALSE)
            Yhat1=pred("rf",cbind(X[,c(NX)],Xbck),X[,NY],cbind(XX1[,c(NX)],Xbck),class=FALSE)
          }
          
          for (k in 1:20){
            Xb<-Xbck[sample(1:N)]
            if (length(sX)>0){
              Yhat0=cbind(Yhat0,pred("rf",cbind(X[,c(sX,NX)],Xb),X[,NY],cbind(XX0[,c(sX,NX)],Xb),class=FALSE))
              Yhat1=cbind(Yhat1,pred("rf",cbind(X[,c(sX,NX)],Xb),X[,NY],cbind(XX1[,c(sX,NX)],Xb),class=FALSE))
            }else{
              Yhat0=cbind(Yhat0,pred("rf",cbind(X[,c(NX)],Xb),X[,NY],cbind(XX0[,c(NX)],Xb),class=FALSE))
              Yhat1=cbind(Yhat1,pred("rf",cbind(X[,c(NX)],Xb),X[,NY],cbind(XX1[,c(NX)],Xb),class=FALSE))
            }       
          }
          Yhat0=apply(Yhat0,1,mean)
          Yhat1=apply(Yhat1,1,mean)
        }
        
        
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))/abs(DD)
        ACC2<-c(ACC2,acc)
        NMSE2<-c(NMSE2, mean(abs(e)))
        cat("Cnt Acc2=",mean(ACC2), "NMSE=",mean(NMSE2), "\n")
        
        ###############
        ############@##### optimal CF
        
        Yhat0=pred("rf",X[,PaY],X[,NY],Xc0[,PaY],class=FALSE)
        Yhat1=pred("rf",X[,PaY],X[,NY],Xc1[,PaY],class=FALSE)
        fs<-PaY
        DX<-NULL
        DY<-NULL
        
        for (i in 1000){
          a=sample(1:NROW(X),2)
          DX<-rbind(X[a[1],fs]-X[a[2],fs])
          DY<-c(DY, X[a[1],NY]-X[a[2],NY])
        }
        
        DYhat=pred("rf",DX,DY,Xc0[,fs]-Xc1[,fs],class=FALSE)
        
        
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        acc2=length(which(((Y0-Y1)*(DYhat))>0))/length(Y0)
        
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))/abs(DD)
        BACC<-c(BACC,acc)
        BACC2<-c(BACC2,acc2)
        BNMSE<-c(BNMSE, mean(abs(e)))
        cat("Best Acc=",mean(BACC),mean(BACC2), "NMSE=", mean(BNMSE), "\n")
        if (acc==0)
          browser()
        
        #########################################
        ############# ALLVAR
        
        
        
        Yhat0=pred("rf",X[,-NY],X[,NY],XX0[,-NY],class=FALSE)
        Yhat1=pred("rf",X[,-NY],X[,NY],XX1[,-NY],class=FALSE)
        
        DX<-NULL
        DY<-NULL
        fs<-setdiff(1:NCOL(X),NY)
        for (i in 1000){
          a=sample(1:NROW(X),2)
          DX<-rbind(X[a[1],fs]-X[a[2],fs])
          DY<-c(DY, X[a[1],NY]-X[a[2],NY])
        }
        
        DYhat=pred("rf",DX,DY,XX0[,fs]-XX1[,fs],class=FALSE)
        
        
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        acc2=length(which(((Y0-Y1)*(DYhat))>0))/length(Y0)
        
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))/abs(DD)
        
        WACC<-c(WACC,acc)
        WACC2<-c(WACC2,acc2)
        WNMSE<-c(WNMSE, mean(abs(e)))
        cat("rr=", rr, " Naive Acc=",mean(WACC), mean(WACC2), "NMSE=", mean(WNMSE), "\n ---\n")
        
        #browser()
        
        
      } # for rr
    }
  } # if ! is.null
}


