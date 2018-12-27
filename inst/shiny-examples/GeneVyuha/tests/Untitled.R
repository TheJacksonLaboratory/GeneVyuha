ui <- fluidPage(
  radioButtons("dist", "Distribution type:",
               c("Normal" = "norm",
                 "Uniform" = "unif",
                 "Log-normal" = "lnorm",
                 "Exponential" = "exp")),
  plotOutput("distPlot")
)

server <- function(input, output) {
  output$distPlot <- renderPlot({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)

    hist(dist(500))
  })
}

fluidRow(
  uiOutput("selectNetworkFile"),
  actionButton("updateSelectedNetwork", "Update")
),

fluidRow(
  textOutput("check")),
fluidRow(
  DTOutput("tb1")),
fluidRow(
  (visNetworkOutput("network1"))),

output$selectNetworkFile <- renderUI({
  networkInFiles <- readRDS("data/storedNetworks.RDS")
  selectInput("selectedNetworkFile", "Try a demo network",
              networkInFiles,
              selected = "None")
})
output$check <- renderText(paste("Notupdated"))
observeEvent(input$updateSelectedNetwork,{
  output$check <- renderText("updated")
  #  if(input$selectedNetworkFile == "None"){

  # }  else {
  ftpo3 <- reactive({
    tmp2 <- read.table(paste0("inputs/",input$selectedNetworkFile, ".tpo"),header = TRUE,
                       stringsAsFactors = FALSE)
    return(tmp2)
  })
  #  }
})



rs2 <- reactive({
  rs2 <-  new("racipeSet")
  network(rs2) <- f_tpo3()
  annotation(rs2) <- input$filenameTopo

  output$tb1 <-renderDT({

    return(f_tpo3())
  }, selection = 'none', editable = FALSE, rownames= FALSE
  )

  output$network1 <- renderVisNetwork({

    plotRSet(rs2(),"network")
  })

  return(rs2)
})

