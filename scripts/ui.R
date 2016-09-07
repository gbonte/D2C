
shiny::shinyUI(fluidPage(
  titlePanel("D2C graphical demo"),
  sidebarPanel(
    fileInput('file1', 'Choose D2C model',
              accept = c(
                '.Rdata'              )
    ),
    textOutput("D2C")),
  column(9, wellPanel(
    plotOutput("graphPlot", height=450)
  )
  ),
  
  column(4,wellPanel(
    sliderInput(inputId = "nNode", label = "nNodes",  value = 15, min=1, max=100),
    sliderInput(inputId = "maxPar", label = "maxPar",  value = 2, min=1, max=10),
    h3(textOutput("nEdges"))))
  ,
  
  column(4,wellPanel(
    sliderInput(inputId = "nSamples", label = "nSamples",  value = 50, min=20, max=1000))
  ),
  column(4,wellPanel(
    sliderInput(inputId = "vertexSize", label = "Vertex Size",  value = 15, min=1, max=100),
    sliderInput(inputId = "arrowSize", label = "Arrow Size",  value = 10, min=1, max=20)
  ))
  ,
  
  column(6, wellPanel(
    actionButton("do", "Assess D2C for this graph  "),
    radioButtons("type", label=" ",
                 choices = list("is.parent" = "is.parent", "is.ancestor" = "is.ancestor", "is.descendant" = "is.descendant", "is.mb" = "is.mb"), 
                 selected = "is.parent"),
    textOutput("text1"),
    textOutput("text2"),
    h3(textOutput("BER")),
    h3(textOutput("AUC")),
    h3(tableOutput('table'))
  ))
  
)
)
