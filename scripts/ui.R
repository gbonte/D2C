shiny::shinyUI(fluidPage(
  titlePanel("igraph plot demo"),
  column(9, wellPanel(
    plotOutput("graphPlot", height=450)
  )
  ),

  column(4, wellPanel(
    numericInput(inputId = "parent", label = "parent",  value = 15, min=1, max=100),
    numericInput(inputId = "son", label = "son",  value = 15, min=1, max=100),
    sliderInput(inputId = "nNode", label = "nNode",  value = 15, min=1, max=100),
    sliderInput(inputId = "nSamples", label = "nSamples",  value = 15, min=10, max=200),
    sliderInput(inputId = "vertexSize", label = "Vertex Size",  value = 15, min=1, max=100),
    sliderInput(inputId = "arrowSize", label = "Arrow Size",  value = 10, min=1, max=20)
  )),
  column(4, wellPanel(
    actionButton("do", "Test"),
    textOutput("text1"),
    textOutput("text2")
  ))

))
