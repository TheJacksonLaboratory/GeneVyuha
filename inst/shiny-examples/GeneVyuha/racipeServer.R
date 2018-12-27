
###########################################
# RACIPE
###########################################



shinyRacipeNetwork <- callModule(shinyLoadNetwork, "shinyRacipeNetwork", stringsAsFactors = FALSE)
observeEvent(input$simulateRacipe, {
  rsRacipe <-new("racipeSet")
  network(rsRacipe) <- shinyRacipeNetwork()
  anneal <- switch(input$sRacipeOption,
                   constantNoise = FALSE,
                   annealing = TRUE,
                   FALSE)

  #      rs <- shinyModelExplorer()
  rsRacipe <- simulateRS(rsRacipe)
  output$racipePca <- renderPlot({
    if(input$simulateRacipe == 0) return()
    plotRSet(rsRacipe, "pca")
  })
  output$racipeHeatmap <- renderPlot({
    if(input$simulateRacipe == 0) return()
    plotRSet(rsRacipe, "exprsHeatmap")
  })

  if(!is.null(rsRacipe)){
    parameterNames <- varLabels(rsRacipe)
    parameters <- params(rsRacipe)
  }
  output$filteredOutputRacipe <- renderUI({
    if(is.null(parameters)) return(NULL)

    selectInput("selectedParameter", "Parameter",
                parameterNames,
                selected = NULL)
  })

  dataSimulation <- normalizeExprs(rsRacipe)
  pca = prcomp(dataSimulation, scale. = F, center = F)


  filtered <- reactive({
    if (is.null(dataSimulation)) {
      return(NULL)
    }
    dataSimulation[((parameters[,input$selectedParameter] >= 0.01*(input$parameterInput[1]*max(parameters[,input$selectedParameter]))) &(parameters[,input$selectedParameter] < 0.01*(input$parameterInput[2]*max(parameters[,input$selectedParameter])))),]

  })


  output$filterSliderRacipe <- renderUI({
    if(is.null(parameters)) return(NULL)
    sliderInput("parameterInput", "Parameter Range", min = 0,
                max = 100, value = c(0,100))
  })

  output$racipePcaFiltered <- renderPlot({
    if(is.null(filtered())) return(NULL)
    pcaData <- scale(filtered(), pca$center, pca$scale) %*% pca$rotation
    pcaData <- data.frame(PC1=pcaData[,1],PC2=pcaData[,2])
    rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
    plotColor <- rf(32)
    binCount <- 40
    plotDensity(pcaData,binCount,plotColor)
  })

  output$racipeHeatmapFiltered <- renderPlot({
    if(is.null(filtered())) return(NULL)
    plotExprsHeatmap(rsRacipe, filtered())
  })

  output$downloadDataRacipe <- downloadHandler(
    filename = function() {
      paste(annotation(rsRacipe), '.RDS', sep='')
    },
    content = function(con) {
     # if(is.null(filtered())) return(NULL)
      if(is.null(rsRacipe)) return(NULL)
      saveRDS(rsRacipe, con)
   #   write.csv(filtered(), con)
    }
  )

})
