
rm(list=ls())
library(devtools)
#install_github("gbonte/D2C")
library(D2C)
library(graph)
library(igraph)
library(gRbase)
library(ROCR)
library(RBGL)
library(bnlearn)




set.seed(0)

type="is.ancestor"
trainD2C<-NULL
knocked<-NULL



additive=FALSE
aBER=NULL
aBER2=NULL
aAUC=NULL
allDAG=NULL
nfeat=40
trainD2C=NULL
for (iter in 0:10000){
    
    print(iter)
    unbalanced=TRUE

    cnt2=1  
    
     
    cnt2=cnt2+1
    set.seed(iter+cnt2)
    maxPar=sample(2:4,1)
    wgt = runif(1,0.8,0.9)
    if (iter %% 10 > 0 || is.null(trainD2C)){
        nNode=cnt2+sample(10:50,1)
        g<-random_dag(1:nNode,maxpar=min(nNode,maxPar),wgt)
        cnt<-0
        
        G<-graph.adjacency(as(g,"matrix"))
        DAGg=as_graphnel(G)
    } else {
        DAGg=allDAG[[which.max(aBER)]]
        G=igraph.from.graphNEL(DAGg)
        
        cat("worst=",which.max(aBER),":",max(aBER),"\n")
        aBER[which.max(aBER)]=-1
        
    }
    allDAG=c(allDAG,DAGg)
    
    if (runif(1)<0.5){
        H = function() return(H_Rn(2)) #function() return(H_sigmoid(1))
    } else {
        if (runif(1)<0.5)
            H = function() return(H_Rn(1))
        else
            H = function() return(H_sigmoid(1))
    }
    
    
    
    const=TRUE

    while (const){
    nSamples=sample(50:500,1)
    
    additive=sample(c(TRUE,FALSE),1)
    weights=c(0.5,1)
    weights[1]=runif(1)
    sdn=runif(1,0.2,0.5)
    DAG = new("DAG.network",
              network=DAGg,H=H,additive=additive,
              weights=weights,sdn=sdn)
    
    observationsDAG <- compute(DAG,N=nSamples)
    const=any(apply(observationsDAG,2,sd)<0.01)
    }
    cat("iter=",iter,"cnt2=",cnt2,"Observed data=",dim(observationsDAG), 
        "nEdges=",length(E(G)),"\n")
    
    trainDAG<-new("simulatedDAG",NDAG=0)
    
    
    trainDAG@NDAG=1
    trainDAG@list.DAGs[[1]]=as_graphnel(G)
    trainDAG@list.observationsDAGs[[1]]=observationsDAG
    
    
    D2C<-new("D2C",sDAG=trainDAG,
             descr=descr,ratioEdges=0.85,
             max.features=nfeat, type=type,goParallel=FALSE,
             verbose=TRUE,npar=1)
    
    N0=length(which((D2C@Y==0)))
    N1=length(which((D2C@Y==1)))
    
    
    print(D2C@origX[1:3,1:3])
    
    
    if (!is.null(trainD2C)){
        Nodes=nodes(DAGg)
        max.edges<-length(edgeList(DAGg))
        
        subset.edges = matrix(unlist(edgeList(DAGg)),ncol=2,byrow = TRUE)
       # subset.edges = unique(rbind(subset.edges,t(replicate(n =max.edges ,
        #                                                     sample(Nodes,size=2,replace = FALSE)))))
        
        
        if (type=="is.parent"){
            subset.edges = unique(rbind(subset.Edges,t(replicate(n =4*max.edges ,sample(Nodes,size=2,replace = FALSE)))))
        } else {
          added=0
          for (n1 in Nodes){
            for (n2 in setdiff(Nodes,n1)){
              if (type=="is.ancestor"){
                if (is.ancestor(DAGg,n1,n2) & (!is.parent(DAGg,n1,n2))){
                  subset.edges = rbind(subset.edges,c(n1,n2))
                  added=added+1
                  cat("+") 
                }
                
              }
              if (type=="is.descendant"){
                if (is.descendant(iDAG2,n1,n2) & (!is.parent(iDAG2,n2,n1)) ){
                  subset.edges = rbind(subset.edges,c(n1,n2)) 
                  added=added+1
                  cat("+")
                }
              }
              if (type=="is.mb"){
                if (is.mb(iDAG2,n1,n2)  )
                  added=added+1
                subset.edges = rbind(subset.edges,c(n1,n2)) 
              }
              
              if (added>sz)  
                break;
            }
          }
        }
        
        Ahat.IAMB<-(amat(iamb(data.frame(observationsDAG),alpha=0.01,max.sx=3)))
        igraph.IAMB<-igraph::graph.adjacency(Ahat.IAMB)
        
        
        Yhat.D2C<-NULL
        Yhat.IAMB<-NULL
        phat.D2C<-NULL
        Ytrue<-NULL
        cat("D2C inferring", NROW(subset.edges), 
            "direct dependencies: please wait \n")
        
        for(jj in 1:NROW(subset.edges)){
            i=subset.edges[jj,1]
            j=subset.edges[jj,2]
            I =as(subset.edges[jj,1],"numeric")
            J =as(subset.edges[jj,2],"numeric") 
            if (length(intersect(c(I,J),knocked))==0){
                pred.D2C.rr =predict(trainD2C,I,J, observationsDAG,rep=4)$prob
                Yhat.D2C<-c(Yhat.D2C,round(pred.D2C.rr))
                phat.D2C<-c(phat.D2C,pred.D2C.rr)
                Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,I,J,type))
                Ytrue<-c(Ytrue,is.what(G,i,j,type)) 
                cat(".")
            }

            i=subset.edges[jj,2]
            j=subset.edges[jj,1]
            I =as(subset.edges[jj,2],"numeric")
            J =as(subset.edges[jj,1],"numeric")
            if (length(intersect(c(I,J),knocked))==0){
                pred.D2C.rr =predict(trainD2C,I,J, observationsDAG,rep=4)$prob
                Yhat.D2C<-c(Yhat.D2C,round(pred.D2C.rr))
                phat.D2C<-c(phat.D2C,pred.D2C.rr)
                Yhat.IAMB<-c(Yhat.IAMB,is.what(igraph.IAMB,I,J,type))
                Ytrue<-c(Ytrue,is.what(G,i,j,type))
                cat(".")
            }

            
        }
        
        
        aBER=c(aBER,BER(Ytrue,Yhat.D2C))
        aAUC=c(aAUC,AUC(Ytrue,phat.D2C))
        aBER2=c(aBER2,BER(Ytrue,Yhat.IAMB))
        cat("\n Assessment accuracy: \n BER D2C=",
            round(mean(aBER),2),"\n")
        cat("AUC D2C=",round(mean(aAUC),2),"\n")
        
        cat("\n  BER IAMB=",
            round(mean(aBER2),2),"\n")
        
    }
    
    if (iter==0){
        allD2C<-D2C
    } else{
        allD2C@origX<-rbind(allD2C@origX,D2C@origX)
        allD2C@Y<-c(allD2C@Y,D2C@Y)
        cat("dim(allD2C@origX)=",dim(allD2C@origX),"NDAG computed=",iter,"\n")
    }
    
    N0=length(which((allD2C@Y==0)))
    N1=length(which((allD2C@Y==1)))

    if (N0>12 & N1>12){
    trainD2C<-makeModel(allD2C,classifier="RF",EErep=4,verbose=FALSE)
    } else {
        trainD2C<-NULL
        }
    cat("D2C learner: \n # descriptors=",NCOL(allD2C@origX), 
        "\n # samples=",NROW(allD2C@origX), "\n # positives=",
        length(which(allD2C@Y==1)) , "\n")
    namefile=paste("./data/sequential",type,"Rdata",sep=".")
    save(file=namefile,list=c("type","allDAG","aBER","aAUC","trainD2C"))
    
}


