## I still keep load these two package for stand alone shiny apps

rm(list=ls())
options(shiny.maxRequestSize=250*1024^2)
library(shiny)
library(graph)
library(igraph)
library(readxl)
library(gRbase)
library(D2C)
library(ROCR)
library(RBGL)
G<-NULL
observationsDAG<-NULL
maxpar<-0
#load(paste("/Users/bontempi/Dropbox/bontempi_office/Rlang/d2c/D2C/data/trainD2C.500","is.mb","RData",sep="."))
#trainD2C.mb<- trainD2C


#load(paste("/Users/bontempi/Dropbox/bontempi_office/Rlang/d2c/D2C/data/trainD2C.1000","is.parent","RData",sep="."))
trainD2C<-NULL



shinyServer(function(input, output) {
  output$GraphTypeUI <- renderUI({
    if (is.null(input$GraphType))
      return()
    
  })
  
  
  
  
  
  
  plotGraph <- function(){
    
    wgt = 0.9
    if (is.null(G) || length(V(G)) != input$nNode || maxpar != input$maxPar){
      g<-random_dag(1:input$nNode,maxpar=min(input$nNode,input$maxPar),wgt)
      
      cnt<-2
      
      while (sum(unlist(lapply(graph::edges(g),length)))<input$nNode & cnt<100){
        g<-random_dag(1:input$nNode,maxpar =min(input$nNode,input$maxPar),wgt)
        cnt<-cnt+1
        
      }
      G<<-graph.adjacency(as(g,"matrix"))
      maxpar<<-input$maxPar
      output$nEdges<-renderText({
        if (!is.null(G))
          paste("nEdges=",length(E(G)))
      })
    }
    
    
    if (input$nSamples != NROW(observationsDAG) || input$nNode != NCOL(observationsDAG)){
      
      if (runif(1)<0.5){
        H = function() return(H_sigmoid(1))
      } else {
        H = function() return(H_Rn(1))
      }
      
      additive=sample(c(TRUE,FALSE),1)
      
      DAG = new("DAG.network",
                network=as_graphnel(G),H=H,additive=additive,sdn=runif(1,0.2,0.5))
      
      observationsDAG <<- compute(DAG,N=input$nSamples)
      
      
      print(dim(observationsDAG))
      
    }
    # Adjust vertex size according user input
    V(G)$size = input$vertexSize
    
    # Adjust arrow size according user input
    E(G)$arrow.size = input$arrowSize/10
    
    
    
    
    # To avoid plot without boundary error
    if(vcount(G) > 0){
      plot(G)
    }
  }
  
  output$graphPlot <- renderPlot({
    # plotGraph()
    suppressWarnings(plotGraph())
  })
  
  
  output$D2C<-renderText({
    if (!is.null(input$file1)){
      L<-load(input$file1$datapath)
      
      paste("# descriptors=",NCOL(trainD2C@X), "\n # samples=",NROW(trainD2C@X), "\n # positives=",length(which(trainD2C@Y==1)) )
    }
  })
  
  observeEvent(input$do,  {
    if (!is.null(input$file1)){
      load(input$file1$datapath)
      
    }
    if (! is.null(G) && ! is.null(trainD2C)){
      
      
      DAG=as_graphnel(G)
      Nodes=nodes(DAG)
      max.edges<-length(edgeList(DAG))
      
      if (type=="is.parent"){
        subset.edges = matrix(unlist(sample(edgeList(DAG),size = max.edges,replace = F)),ncol=2,byrow = TRUE)
        subset.edges = unique(rbind(subset.edges,t(replicate(n =3*max.edges ,sample(Nodes,size=2,replace = FALSE)))))
      } else {
        subset.edges = unique(t(replicate(n =4*max.edges ,sample(Nodes,size=2,replace = FALSE))))
      }
      
      Yhat.D2C<-NULL
      phat.D2C<-NULL
      Ytrue<-NULL
      for(jj in 1:NROW(subset.edges)){
        i=subset.edges[jj,1]
        j=subset.edges[jj,2]
        I =as(subset.edges[jj,1],"numeric");
        J =as(subset.edges[jj,2],"numeric") ;
        pred.D2C = predict(trainD2C,I,J, observationsDAG)
        
        Yhat.D2C<-c(Yhat.D2C,as.numeric(pred.D2C$response)  -1)
        
        phat.D2C<-c(phat.D2C,pred.D2C$prob[1,"1"])
        Ytrue<-c(Ytrue,is.what(G,i,j,input$type)) ##graphTRUE[subset.edges[jj,1],subset.edges[jj,2]])
        
        cat(".")
        
      }
      output$BER<-renderText({
        paste("BER=",round(BER(Ytrue,Yhat.D2C),2))
      })
      
      output$AUC<-renderText({
        paste("AUC=",round(AUC(Ytrue,phat.D2C),2))
      })
      print(table(Ytrue,round(Yhat.D2C)))
      output$table <- renderTable({
        
        
        A=table(Ytrue,round(Yhat.D2C))
        rownames(A)=c("N","P")
        colnames(A)=c("N'","P'")
        A
        #paste("BER",BER.D2C)
      })
    }
    
    
  })
  
  ## Output of Adjacency Matrix Panel
  calculateAdjMatrix <- function(){
    g=realTimeGraph()
    if(is.weighted(g)){
      adjmat <- as(as_adjacency_matrix(g, attr='weight'),"matrix")
    }else{
      adjmat <- as(as_adjacency_matrix(g),"matrix")
    }
    return(adjmat)
  }
  output$AdjMatrix <- renderTable(calculateAdjMatrix())
  
  # output$AdjMatrix <- renderTable(as(as_adjacency_matrix(realTimeGraph(), attr='weight'),"matrix"))
  
  ## Output of Centrality Panel
})
