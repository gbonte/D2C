## I still keep load these two package for stand alone shiny apps
library(shiny)
library(graph)
library(igraph)
library(readxl)
library(gRbase)
library(D2C)

G<-NULL
observationsDAG<-NULL

load(paste("/Users/gbonte/Dropbox/bontempi_office/Rlang/d2c/D2C/data/trainD2C.500","is.mb","RData",sep="."))
trainD2C.mb<- trainD2C


load(paste("/Users/gbonte/Dropbox/bontempi_office/Rlang/d2c/D2C/data/trainD2C.200","is.parent","RData",sep="."))


shinyServer(function(input, output) {
  output$GraphTypeUI <- renderUI({
    if (is.null(input$GraphType))
      return()

  })






  plotGraph <- function(){

    wgt = runif(n = 1,min = 0.65,max = 0.95)
    if (is.null(G) || length(V(G)) != input$nNode){
      g<-random_dag(1:input$nNode,maxpar=4,wgt)
      G<<-graph.adjacency(as(g,"matrix"))

    }


    if (input$nSamples != NROW(observationsDAG)){

      H = function() return(H_sigmoid(1))

      DAG = new("DAG.network",
                network=as_graphnel(G),H=H,additive=TRUE,sdn=0.2)
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

  observeEvent(input$do, {

    if (input$parent == input$son){
      output$text1 <- renderText({
        "Son equal to parent"
      })
    } else {
      if (input$parent > input$nNode ||  input$son > input$nNode){
        output$text1 <- renderText({
          "Too large values"
        })
      } else {
        p<-predict(trainD2C,input$parent,input$son, observationsDAG)
        p2<-predict(trainD2C,input$son,input$parent, observationsDAG)

        p.mb<-predict(trainD2C.mb,input$parent,input$son, observationsDAG)
        p2.mb<-predict(trainD2C.mb,input$son,input$parent, observationsDAG)

        print(p$prob[1,"1"])
        output$text1 <- renderText({
          paste("->",p$prob[1,"1"]," <-", p2$prob[1,"1"])
        })

        output$text2 <- renderText({
          paste("->",p.mb$prob[1,"1"]," <-", p2.mb$prob[1,"1"])
        })

      }
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
