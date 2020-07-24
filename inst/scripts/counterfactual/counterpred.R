## D2CF counterfactual
## Script associated to the article 
## "Beyond uncounfoundness in predicting counterfactuals: a machine learning approach"
## by G. Bontempi

rm(list=ls())
library(randomForest)
library(graph)
library(devtools)
#install_github("gbonte/D2C")
#install_github("gbonte/gbcode")
library(D2C)
library(gbcode)
load("data/sequential.Rdata")
trainD2CP<-trainD2C
load("./data/sequential.is.descendant.Rdata")
trainD2CD<-trainD2C
load("./data/sequential.is.ancestor.Rdata")

verbose=FALSE

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
  
  set.seed(r)
  
  n=sample(10:20,1)
  # number of variables
  
  wgt = runif(n = 1,min = 0.85,max = 1)
  
  rDAG<-gRbase::random_dag(1:n,maxpar = 3,wgt)
  nodes(rDAG)<-as.character(1:n)
  
  Edges= gRbase::edgeList(rDAG) 
  Nodes<-nodes(rDAG)
  irDAG=igraph::graph.adjacency(as(rDAG,"matrix"))
  
  bestSP=0
  NX=NULL
  listN=NULL
  for (nn in Nodes){
    for (jj in setdiff(Nodes,nn)){
      SP=as.numeric(igraph::shortest.paths(irDAG, nn, jj, "out"))
      if ((!is.infinite(SP)) & SP >0 ){
        
        NX=nn
        NY=jj
        nascX=0
        for (nnn in Nodes){
          if (is.descendant(irDAG,nnn,NX)  ){
            nascX<-nascX+1     
          }}
        
        if (  SP>2 & nascX>3 )
          listN=c(listN,list(c(NX,NY,SP,nascX)))
        
        
      }
    }
  } ## for nn
  
  ## listN: list of tuples (NX,NY) where NX is the traitment and NY is the outcome
  
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
      
      if (verbose)
      cat("r=", r, "rr=",rr, "NX=",NX,"NY=", NY,"LP=",listN[[rr]][3],
          "NascX=",listN[[rr]][4], "PaY=",PaY, "ChY=",ChY,"DescX=",DescX,"\n")
      NY<-as.numeric(NY) 
      NX=as.numeric(NX)
      H1 = function() return(H_Rn(1))
      H2 = function() return(H_Rn(2))
      Hs = function() return(H_sigmoid(1))
      Hk = function() return(H_kernel())
      
      cnt=0
      L=0
      DY=0
      
      ## Loop to create a counterfactual dataset such that the treatment variation has
      ## sufficient impact on outcome 
      
      while (L<100 | max(abs(DY))>4){
        ## max(abs(DY)): it avoids ill conditioned configurations with excessive outcome variation 
        maxN=sample(50:1000,1)
        set.seed(r+cnt)
        cnt=cnt+1
        HH=sample(c(Hs,H1,H2),1)[[1]]
        DAG = new("DAG.network",network=rDAG,
                  additive=sample(c(TRUE,FALSE),1),
                  sdn=c(0.1,0.2),exosdn=runif(1,0.15,2),
                  weights=c(0.5,2) ,H = HH)
        
        X=compute(DAG,maxN) 
        N=NROW(X)
        
        if (N>(maxN/3)){ ## check on the number of simulated samples
          if ( all(apply(X,2,sd)>0.01)) { ## check to avoid constant variables
            set.seed(as.numeric(Sys.time())+cnt)
            DX=runif(1,0.5,2)
            
            ## control
            delta0=sample(X[,NX],N,rep=TRUE)
            Xc0=counterfact(DAG,X,delta=delta0, knocked=NX)
            Y0=Xc0[,NY]
            
            ## case
            delta1=X[,NX]+runif(N,DX,2*DX)*sample(c(-1,1),N,rep=TRUE)
            Xc1=counterfact(DAG,X,delta=delta1 ,knocked=NX)
            Y1=Xc1[,NY]
            
            DY=Y1-Y0
            ## counterfactual variability of the outcome
            
            L=length(unique(round(DY,2)))
            ## measure of variability of the outcome
            
            if (L>100 & verbose)
              cat("L=",L,"DX=",DX,"summary(DY)=",summary(DY),"\n")
            
          }
          
        }
        
        if (cnt>10){
          ## too many attempts: let's change  the DAG
          break
        }
      }
      
      
      
      if (L>100 & max(abs(DY))<4){
        X=scale(X)
        Xc0=scale(Xc0,center=attr(X,'scaled:center'),scale=attr(X,'scaled:scale'))
        Xc1=scale(Xc1,center=attr(X,'scaled:center'),scale=attr(X,'scaled:scale'))
        Y0=Xc0[,NY]
        Y1=Xc1[,NY]
        
        
        ancX=numeric(n)
        parY=numeric(n)
        descX=numeric(n)
        
        
        ## prediction by D2C
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
        
        ## D2C prediction: set of variables parent of Y and no descendants of X
        px=which.max(ancX)
        if (max(ancX)>0.5)
          px=sort(ancX,decr=TRUE,index=TRUE)$ix[1:2]
        ## D2C prediction: set of variables parent of X 
        if (verbose)
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
        ############### D2CF counterfactual
        
        if (length(sX)>0){
          Yhat0=pred("rf",X[,c(sX,NX)],X[,NY],XX0[,c(sX,NX)],class=FALSE)
          Yhat1=pred("rf",X[,c(sX,NX)],X[,NY],XX1[,c(sX,NX)],class=FALSE)
        }
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))
        ACC<-c(ACC,acc)
        NMSE<-c(NMSE, mean(abs(e)))
        cat("D2Counterf Accuracy=",mean(ACC), "ABSEerr=",mean(NMSE), "\n")
        
        ############### D2CF MC
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
        e=(DD -(Yhat0-Yhat1))
        ACC2<-c(ACC2,acc)
        NMSE2<-c(NMSE2, mean(abs(e)))
        cat("D2Counterf (bis) Accuracy =",mean(ACC2), "ABSErr=",mean(NMSE2), "\n")
        
        ###############
        ############@##### ORACLE CF
        
        Yhat0=pred("rf",X[,PaY],X[,NY],Xc0[,PaY],class=FALSE)
        Yhat1=pred("rf",X[,PaY],X[,NY],Xc1[,PaY],class=FALSE)
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))
        BACC<-c(BACC,acc)
        
        BNMSE<-c(BNMSE, mean(abs(e)))
        cat("ORACLE Accuracy=",mean(BACC), "ABSErr=", mean(BNMSE), "\n")
        
        #########################################
        ############### Naive single input (use only X)
        
        Yhat0=pred("rf",X[,NX],X[,NY],XX0[,NX],class=FALSE)
        Yhat1=pred("rf",X[,NX],X[,NY],XX1[,NX],class=FALSE)
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))
        ACC0<-c(ACC0,acc)
        NMSE0<-c(NMSE0, mean(abs(e)) )
        cat("Naive (single input) Accuracy =",mean(ACC0), 
            "ABSErr=",mean(NMSE0), "\n")
        
        ############# Naive ALLVAR
        Yhat0=pred("rf",X[,-NY],X[,NY],XX0[,-NY],class=FALSE)
        Yhat1=pred("rf",X[,-NY],X[,NY],XX1[,-NY],class=FALSE)
        acc=length(which(((Y0-Y1)*(Yhat0-Yhat1))>0))/length(Y0)
        
        DD=  (Y0-Y1)
        e=(DD -(Yhat0-Yhat1))
        
        WACC<-c(WACC,acc)
        WNMSE<-c(WNMSE, mean(abs(e)))
        cat("Naive (all inputs) Accuracy=",mean(WACC),  
            "ABSErr=", mean(WNMSE), "\n ---\n")
        
        save(file="./data/res.counter.Rdata",list=c("ACC0","NMSE0","WACC","WNMSE",
                                                    "ACC2","NMSE2","ACC","NMSE",
                                                    "BACC","BNMSE"))
        
        
      } # for rr
    }
  } # if length(listN)
} ## for r


